function [TotFailedBanks,Norm_TotFailedBanks,TotCapitalLoss,Norm_TotCapitalLoss,banks_CDM_results,shockterm] =...
    cascadealgorithm(banks,n_banks,shockedbanks,ActiveBanks,FailedBanks,IB_loan_matrix,IB_shockprop_matrix,Tot_init_capital,shockpars,dispoutput)

banks_CDM_loss     = zeros(n_banks,1);
banks_CDM_failed   = zeros(n_banks,1);
banks_CDM_failtime = zeros(n_banks,1);

shockmagnitude  = shockpars(1);

t = 1;

for i = ActiveBanks
    banks(i).shock(t) = 0;
    banks(i).failtime = 0;
    banks(i).failed   = 0;  
end

banks(shockedbanks).shock(t) = banks(shockedbanks).assets.externalassets(t)*shockmagnitude;
banks(shockedbanks).assets.externalassets(t+1) = banks(shockedbanks).assets.externalassets(t) - banks(shockedbanks).shock(t);

% (II) Continue propagation until all banks have failed

while numel(FailedBanks) ~= n_banks  && t <= 1000
    
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
            banks(i).failtime  = t;
            banks(i).failed    = 1;
                    
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
    
    shockedbanks_update = unique(shockedbanks_update); % Control for the fact that shocked banks can have creditors in common
    
    if absorbcount == numel(shockedbanks) % Break loop prematurely if shock is absorbed by all banks in current propagation round
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
    shockedbanks = shockedbanks_update; % Update set of shocked banks to pass onto next iteration round
    t=t+1;
end

for i = 1:n_banks
    lossvec(i,:) = banks(i).loss;
    
    banks_CDM_loss(i)      = sum(banks(i).loss);
    banks_CDM_failed(i)   = banks(i).failed;
    banks_CDM_failtime(i) = banks(i).failtime;
    
end

banks_CDM_results = [banks_CDM_loss banks_CDM_failed banks_CDM_failtime];

% Dependant variables
TotFailedBanks      = sum(NumFailedBanks);
Norm_TotFailedBanks = TotFailedBanks/n_banks; % Normalisation: Divide Failed banks by total number of banks

TotCapitalLoss      = sum(sum(lossvec)); 
Norm_TotCapitalLoss = TotCapitalLoss/Tot_init_capital; % Normalisation: Divide capital loss by initial system capital






end