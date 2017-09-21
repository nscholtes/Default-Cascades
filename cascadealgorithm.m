function [banks,TotFailedBanks,Norm_TotFailedBanks,TotCapitalLoss,Norm_TotCapitalLoss,banks_CDM_results,shockterm,initshock_magnitude,assetprice] =...
    cascadealgorithm(banks,n_banks,shockedbanks,ActiveBanks,FailedBanks,probability_matrix,...
    IB_loan_matrix,IB_shockprop_matrix,Tot_init_capital,Tot_init_ext_assets,shockpars,dispoutput,liquidityrisk,tstop)

banks_CDM_loss     = zeros(n_banks,1);
banks_CDM_failed   = zeros(n_banks,1);
banks_CDM_failtime = zeros(n_banks,1);

assetprice = 1;

shockmagnitude  = shockpars(1);
omega           = shockpars(3);

system_ext_assets = Tot_init_ext_assets;

t = 1;

for i = ActiveBanks
    banks(i).shock(t) = 0;
    banks(i).failtime = 0;
    banks(i).failed   = 0;  
end

% Loop to reduce external assets of selected initially shocked bank(s)
for i = shockedbanks

    banks(i).shock(t) = banks(i).assets.externalassets(t)*shockmagnitude;
    banks(i).assets.externalassets(t+1) = banks(i).assets.externalassets(t) - banks(i).shock(t);

    initshock_magnitude = sum(banks(i).shock(t));
end

% (II) Continue propagation until all banks have failed (with break loop conditions in lines 198 & 204)

while numel(FailedBanks) ~= n_banks  && t <= tstop
    
    failurecount = 0; % Count number of failures in current propagation step
    absorbcount  = 0; % Count number of times shock is absorbed by bank capital (no further propagation in this case)
    
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
        
        banks(i).num_lenders(t+1)   = banks(i).num_lenders(t);
        banks(i).num_borrowers(t+1) = banks(i).num_borrowers(t);
            
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
         
        % CASE I: Initially shocked bank's capital sufficient to absorb loss - No further contagion from current shocked bank
        if banks(i).shock(t)  < banks(i).liabilities.capital(t)
            
            absorbcount = absorbcount+1;
            
            banks(i).status(t)    = 1;
            banks(i).firesales(t) = 0;
            
            if t > 1
                banks(i).assets.externalassets(t+1) = banks(i).assets.externalassets(t);
            end
            
            banks(i).liabilities.capital(t+1)   = banks(i).liabilities.capital(t) - banks(i).shock(t); % Reduce capital by shock amount
            
            banks(i).loss(t) = banks(i).shock(t); 
            
            if dispoutput == 1
                fprintf(1,'Bank %d is able to absorb the shock in period %d with a capital loss of %.3f\n',i,t,banks(i).shock(t));
            end 
            
        % CASE II: Insufficient capital to absorb loss. Bank i fails and propagates to connected creditors 
        elseif banks(i).shock(t) >= banks(i).liabilities.capital(t) 
            
            failurecount = failurecount+1;
 
            banks(i).status(t) = 0;
            banks(i).failtime  = t;
            banks(i).failed    = 1;
            
            % Firesale mechanism: Sell off all external assets
            banks(i).assets.externalassets(t+1) = 0;
            banks(i).firesales(t) = banks(i).assets.externalassets(t); 
                    
            shockedbanks_update = [shockedbanks_update banks(i).lender_ids']; % After i's failure, update set of shocked banks = i's creditors
            
            banks(i).liabilities.capital(t+1) = 0;
            
            banks(i).loss(t) = banks(i).liabilities.capital(t);
            
            banks(i).residualshock(t) = banks(i).shock(t) - banks(i).liabilities.capital(t);
            
            Lvec_P = sprintf(' %.3g ',banks(i).lender_ids);
            
            if dispoutput == 1
                if banks(i).num_lenders(t) ~= 0
                    fprintf(1,'Bank %d fails in period %d and firesales %.3f worth of external assets at current market price %.3f\n',...
                        i,t,banks(i).firesales(t),assetprice(t));
                    
                    fprintf(1,'--> Transmits %.3f to %d creditors: [%s]\n',banks(i).residualshock(t),banks(i).num_lenders(t+1),Lvec_P)
                else
                    fprintf(1,'--> Bank %d fails in period %d but has no creditors to transmit shock to\n',i,t);
                end
            end
            
            % Case II-1: Shock absorbed by interbank liabilities               
            if banks(i).residualshock(t)  < banks(i).liabilities.IB_Tot_Borrowing(t)
               
               banks(i).shocktransmit = ones(banks(i).num_lenders(t+1),1)*(1/banks(i).num_lenders(t+1)).*banks(i).residualshock(t);
                                                                                    
            % Case II-2: Residual shock larger than capital AND interbank liabilities --> Customer deposits as ultimate sink
            else         
                banks(i).shocktransmit             = IB_loan_matrix(banks(i).lender_ids,i,t);
                banks(i).residualshock2(t)         = banks(i).residualshock(t) - banks(i).liabilities.IB_Tot_Borrowing(t);
                banks(i).liabilities.deposits(t+1) = 0;           
            end
                        
            % Update loan matrix to account for shock transmission between connected banks
            IB_loan_matrix(banks(i).lender_ids,i,t+1)    = IB_loan_matrix(banks(i).lender_ids,i,t) - banks(i).shocktransmit;            
            IB_shockprop_matrix(banks(i).lender_ids,i,t) = banks(i).shocktransmit;
            
            % Remove failed banks from lender and borrower sets
            FailedBanks    = [FailedBanks i];
            
            ActiveBanks(ActiveBanks == i) = [];

            for j = ActiveBanks
            	banks(j).lender_ids(banks(j).lender_ids == i)   = [];
                banks(j).num_lenders(t+1) = numel(banks(j).lender_ids);
                
                banks(j).borrower_ids(banks(j).borrower_ids == i) = [];
                banks(j).num_borrowers(t+1) = numel(banks(j).borrower_ids);       
            end          
        end
    end
        
    for i = 1:numel(shockedbanks)
        FSvec(i) = banks(shockedbanks(i)).firesales(t);
    end
    
    totalfiresales(t) = sum(FSvec);
    
    clearvars FSvec
    
    % Inverse demand function based on total external asset firesales
    
    if totalfiresales(t) ~=0
        assetprice(t+1) = assetprice(t)*exp(-omega*(totalfiresales(t)/system_ext_assets(t)));
    end
    
    system_ext_assets(t+1) = system_ext_assets(t) - totalfiresales(t);
    
    %assetprice
    
    if dispoutput == 1 && totalfiresales(t)~=0
        disp('-------------------- External asset dynamics -------------------- ')
        fprintf(1,'%.3f worth of external assets were sold in firesale in period %d\n',totalfiresales(t),t);
        fprintf(1,'--> Asset price change: %.3f -> %.3f\n',assetprice(t),assetprice(t+1));
        disp('------------------------------------------------------------------')
    elseif dispoutput == 1 && totalfiresales(t) == 0 
        disp('-------------------- External asset dynamics -------------------- ')
        fprintf(1,'No firesales in period %d\n',t)
        fprintf(1,'--> Asset price does not change from its previous value of %.3f\n',assetprice(t))
        disp('------------------------------------------------------------------')
    end
 
    for i = shockedbanks
        if liquidityrisk == 1 && banks(i).status(t) == 0 % COND: Liquidity risk enabled and shocked bank fails in current period
            banks(i).shocktransmit = banks(i).shocktransmit+...
            ones(banks(i).num_lenders(t+1),1)*(1/banks(i).num_lenders(t+1)).*(assetprice(t) - assetprice(t+1))*banks(i).firesales(t);              
        end
    end
    
    shockedbanks_update = unique(shockedbanks_update); % Control for the fact that shocked banks can have creditors in common
    
    for j = FailedBanks
        shockedbanks_update(shockedbanks_update == j) = [];
    end
    
    if absorbcount == numel(shockedbanks) % Break loop condition 1: If shock is absorbed by all banks in current propagation round
        if dispoutput == 1
            disp('                       o----------------------------------o')
            fprintf('                       | End of shock cascade in period %d |\n',t);
            disp('                       o----------------------------------o')
        end
        shockterm  = t;
        assetprice = assetprice(1:shockterm); 
        break
    end
    
    if isempty(shockedbanks_update) % Break loop  condition 2: Current set of shocked banks do not have creditors for onward propagation
        if dispoutput == 1
        disp('                       o-------------------------------------------------------o')
            fprintf('                       | No more creditors for further propagation in period %d |\n',t);
        disp('                       o-------------------------------------------------------o')       
        end
        shockterm  = t;
        assetprice = assetprice(1:shockterm); 
        break
    end
    
    NumFailedBanks(t) = failurecount;
    
    nonshockedbanks   = setdiff(ActiveBanks,shockedbanks);
    
    % Update values for banks that weren't shocked in the current period
    for i = nonshockedbanks
        banks(i).liabilities.capital(t+1)   = banks(i).liabilities.capital(t);
        banks(i).assets.externalassets(t+1) = banks(i).assets.externalassets(t);
        banks(i).status(t) = 1;
        banks(i).loss(t)   = 0;
    end
    shockedbanks = shockedbanks_update; % Update set of shocked banks to pass onto next iteration round
    t=t+1;
end

for i = 1:n_banks
    lossvec(i,:) = banks(i).loss;
    
    banks_CDM_loss(i)     = sum(lossvec(i,:));
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