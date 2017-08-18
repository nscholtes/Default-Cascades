function [banks,CascadeOutput] = cascadingdefaults(banks,n_banks,ActiveBanks,FailedBanks,adjacency_matrix,probability_matrix,ratios,shockpars,shockspecs)

%--------------------------------------------------------------------------
%% Initialisation and preallocation
%--------------------------------------------------------------------------

theta = ratios(1);
gamma = ratios(2);

shockmagnitude  = shockpars(1);
numshockedbanks = shockpars(2);

shocktype       = shockspecs{1};
targetcriterion = shockspecs{2};

dispoutput = 0; % Show output on propagation dynamics (for debugging and future modifications)

IB_loan_matrix      = zeros(n_banks,n_banks);
IB_shockprop_matrix = zeros(n_banks,n_banks);

%--------------------------------------------------------------------------
%% Calibrating heterogenous bank balance sheets
%--------------------------------------------------------------------------

for i = ActiveBanks    
    
    banks(i).status = 1;
    
    % (I) Counterparty information from adjacency matrix 
    borrower_vec = adjacency_matrix(i,:);
    lender_vec   = adjacency_matrix(:,i);
    
    banks(i).borrower_ids = find(borrower_vec);
    banks(i).lender_ids   = find(lender_vec);
    
    banks(i).num_borrowers = numel(banks(i).borrower_ids);
    banks(i).num_lenders   = numel(banks(i).lender_ids);
    
    % (II) Balance sheet composition from bank size distribution and calibrated parameters
    
    banks(i).assets.externalassets = theta*banks(i).assets.total;
    banks(i).assets.IB_Tot_Lending = (1-theta)*banks(i).assets.total;
    
    banks(i).borrower_peculiarity = probability_matrix(i,banks(i).borrower_ids)';
    
    % (III) Interbank lending from (I) & (II)

    if banks(i).num_borrowers ~= 0
        banks(i).assets.IB_Bil_Lending = (banks(i).assets.IB_Tot_Lending*banks(i).borrower_peculiarity)./(sum(banks(i).borrower_peculiarity));
        
        IB_loan_matrix(i,banks(i).borrower_ids) = banks(i).assets.IB_Bil_Lending;
        banks(i).IB_IsLender = 1;
    else
        banks(i).assets.IB_Tot_Lending = 0;
        banks(i).assets.IB_Bil_Lending = 0;
        banks(i).IB_IsLender    = 0;
        
        banks(i).assets.externalassets = banks(i).assets.total;
    end        
end

%  (IV) Bank borrower information based on (III)

for i = ActiveBanks    
    if banks(i).num_lenders ~=0
        banks(i).liabilities.IB_Bil_Borrowing = IB_loan_matrix(banks(i).lender_ids,i);
        banks(i).liabilities.IB_Tot_Borrowing = sum(banks(i).liabilities.IB_Bil_Borrowing);
        banks(i).IB_IsBorrower = 1;
    else
        banks(i).liabilities.IB_Bil_Borrowing = 0;
        banks(i).liabilities.IB_Tot_Borrowing = 0;
        banks(i).IB_IsBorrower = 0;
    end   
    
    banks(i).IB_IsIntermediary = 0;
end

% (V) Classifying banks who borrow AND lend as intermediaries

for i = ActiveBanks
    if (banks(i).IB_IsLender == 1 && banks(i).IB_IsBorrower == 1)
        banks(i).IB_IsIntermediary = 1;
    end
end

% (VI) Determining bank equity and deposits

for i = ActiveBanks
    banks(i).liabilities.total    = banks(i).assets.total;
    banks(i).liabilities.capital  = gamma*banks(i).liabilities.total;
    banks(i).liabilities.deposits = banks(i).assets.total - (banks(i).liabilities.capital+banks(i).liabilities.IB_Tot_Borrowing);   
end

% Collecting all borrower, lender and intermediary information

for i = ActiveBanks    
    Borrower_vec(i)     = banks(i).IB_IsBorrower;
    Lender_vec(i)       = banks(i).IB_IsLender;
    Intermediary_vec(i) = banks(i).IB_IsIntermediary;    
end

All_Borrowers      = find(Borrower_vec);
All_Lenders        = find(Lender_vec);
All_Intermediaries = find(Intermediary_vec);

% Reorganising bank structure for clarity

permutationvector  = [1 9 3 7 5 4 6 8 10 11 2];
tempstore = banks;

clearvars banks
banks = orderfields(tempstore,permutationvector);

%--------------------------------------------------------------------------
%% Propagation dynamics
%--------------------------------------------------------------------------

%load('test.mat')

% (I) Initial shock to external assets

if strcmp(shocktype,'random')
    shockedbanks = randsample(All_Borrowers,numshockedbanks);
end

shockedbanks = 1;

initial_shockedbank = shockedbanks;

t=1;

for i = ActiveBanks
    banks(i).shock(t) =0;
end

banks(shockedbanks).shock(t) = banks(shockedbanks).assets.externalassets(t)*shockmagnitude;
banks(shockedbanks).assets.externalassets(t+1) = banks(shockedbanks).assets.externalassets(t) - banks(shockedbanks).shock(t);

while numel(FailedBanks) ~= n_banks 
    
    failurecount = 0; % Count number of failures in current propagation step
    absorbcount  = 0; 
    
    if dispoutput == 1
        disp('---------------------------------------------------------------------------------------------')
        fprintf(1,'The current iteration step is %d\n',t);
        disp('---------------------------------------------------------------------------------------------')
    end
    
    IB_loan_matrix(:,:,t+1) = IB_loan_matrix(:,:,t);
    
    for i = ActiveBanks
        banks(i).liabilities.IB_Bil_Lending    = nonzeros(IB_loan_matrix(i,:,t));
        banks(i).liabilities.IB_Tot_Lending(t) = sum(banks(i).liabilities.IB_Bil_Lending);
        
        banks(i).liabilities.IB_Bil_Borrowing      = nonzeros(IB_loan_matrix(:,i,t));
        banks(i).liabilities.IB_Tot_Borrowing(t) = sum(banks(i).liabilities.IB_Bil_Borrowing);
        
        banks(i).status(t) = 1;
        banks(i).loss(t)   = 0;
        
    end
    
    for i = FailedBanks
        banks(i).status(t) = 0;
        banks(i).loss(t)   = 0;
    end

    shockedbanks_update = [];
    
    for i = shockedbanks
         if t > 1
            banks(i).shock(t) = sum(IB_shockprop_matrix(i,:,t-1));
         end
         
        % Case I: Initially shocked bank's capital sufficient to absorb loss - No further contagion from current shocked bank
        if banks(i).shock(t)  < banks(i).liabilities.capital(t)
            
            absorbcount = absorbcount+1;
            
            banks(i).status(t) = 1;
        
            banks(i).liabilities.capital(t+1) = banks(i).liabilities.capital(t) - banks(i).shock(t); % Reduce capital by shock amount
            
            banks(i).loss(t) = banks(i).shock(t); 
            
            if dispoutput == 1
                fprintf(1,'--> Bank %d is able to absorb the shock in period %d with a capital loss of %.3f\n',i,t,banks(i).shock(t));
            end 
            
        % Case II: Insufficient capital to absorb loss. Bank fails and propagates to connected creditors 
        elseif banks(i).shock(t) >= banks(i).liabilities.capital(t) 
            
            failurecount = failurecount+1;
 
            banks(i).status(t) = 0;
                    
            shockedbanks_update = [shockedbanks_update banks(i).lender_ids']; % After i's failure, update set of shocked banks = i's creditors
            
            banks(i).liabilities.capital(t+1) = 0;
            
            banks(i).loss(t) = banks(i).liabilities.capital(t);
            
            banks(i).residualshock(t) = banks(i).shock(t) - banks(i).liabilities.capital(t);
            
            Lvec_P = sprintf(' %.3g ',banks(i).lender_ids);
            
            if dispoutput == 1
                if banks(i).num_lenders ~= 0
                    fprintf(1,'--> Bank %d fails in period %d and trasmits %.3f to %d creditors: [%s]\n',...
                        i,t,banks(i).residualshock(t),banks(i).num_lenders(t),Lvec_P)
                else
                    fprintf(1,'--> Bank %d fails in period %d but has no creditors to transmit shock to\n',i,t);
                end
            end
            
            % Case II-1: Shock absorbed by interbank liabilities
            if banks(i).residualshock(t)  < banks(i).liabilities.IB_Tot_Borrowing(t)
                
               banks(i).shocktransmit = ((1-probability_matrix(banks(i).lender_ids,i))./...
                   sum(1-probability_matrix(banks(i).lender_ids,i)))*banks(i).residualshock(t);
                                           
            % Case II-2: Residual shock larger than capital AND interbank liabilities --> Customer deposits as ultimate sink
            else
                banks(i).shocktransmit             = IB_loan_matrix(banks(i).lender_ids,i,t);
                banks(i).residualshock2(t)         = banks(i).residualshock(t) - banks(i).liabilities.IB_Tot_Borrowing(t);
                banks(i).liabilities.deposits(t+1) = 0;           
            end
            
            % Update loan matrix to account for shock transmission between banks
            IB_loan_matrix(banks(i).lender_ids,i,t+1)    = IB_loan_matrix(banks(i).lender_ids,i,t) - banks(i).shocktransmit;              
            IB_shockprop_matrix(banks(i).lender_ids,i,t) = banks(i).shocktransmit;
            
            % Remove failed banks from lender and borrower sets
            for j = ActiveBanks
            	banks(j).lender_ids(banks(j).lender_ids == i)   = [];
                banks(j).num_lenders(t+1) = numel(banks(j).lender_ids);
                   
                banks(j).borrower_ids(banks(j).borrower_ids == i) = [];
                banks(j).num_borrowers(t+1) = numel(banks(j).borrower_ids);       
            end

            FailedBanks    = [FailedBanks i]; 
            ActiveBanks(ActiveBanks == i) = [];
               
            for j = numel(FailedBanks)
                shockedbanks_update(shockedbanks_update == FailedBanks(j)) = [];
            end            
        end
    end
    
    if absorbcount == numel(shockedbanks)
        if dispoutput == 1
            disp('                       o----------------------------------o')
            fprintf('                       | End of shock cascade in period %d |\n',t);
            disp('                       o----------------------------------o')
        end
        shockterm = t;
        break
     end
    
    NumFailedBanks(t) = failurecount;
    
    nonshockedbanks   = setdiff(ActiveBanks,shockedbanks);
    
    % Update values for banks that weren't shocked in the current period
    for i = nonshockedbanks
        banks(i).liabilities.capital(t+1) = banks(i).liabilities.capital(t);
        banks(i).status(t) = 1;
        banks(i).loss(t)   = 0;
    end
    
    shockedbanks = shockedbanks_update;
    t=t+1;
end

for i = 1:n_banks
    lossvec(i,:) = banks(i).loss;
end

TotFailedBanks = sum(NumFailedBanks);
TotCapitalLoss = sum(sum(lossvec)); 

CascadeOutput = [TotFailedBanks, TotCapitalLoss, shockterm];
%--------------------------------------------------------------------------
% Plot cascading default dynamics
%--------------------------------------------------------------------------
% 
%  t  = 0:0.1:10;
%  y1 = sin(t);
%  y2 = sin(t + 0.1) + 0.5;
%  y3 = cos(t).^2;
%  y4 = cos(t + 0.2).^2 - 0.2;
% 
%  figure
%  hold on
%  fill([t fliplr(t)],[y2 fliplr(y1)] ,'r')
%  fill([t fliplr(t)],[y4 fliplr(y3)] ,'b')
%  grid
%  alpha(0.15)

end