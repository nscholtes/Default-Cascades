function [banks,TotFailedBanks,Norm_TotFailedBanks,TotCapitalLoss,Norm_TotCapitalLoss,TotDepositLoss,Norm_TotDepositLoss,...
    banks_CDM_results,shockterm,tot_initshock_magnitude,assetprice] =...
    cascadealgorithm(banks,n_banks,shockedbanks,ActiveBanks,FailedBanks,probability_matrix,...
    IB_loan_matrix,IB_shockprop_matrix,Tot_init_capital,Tot_init_ext_assets,Tot_init_deposits,shockpars,dispoutput,liquidityrisk,tstop)

banks_CDM_loss     = zeros(n_banks,1);
banks_CDM_failed   = zeros(n_banks,1);
banks_CDM_failtime = zeros(n_banks,1);

assetprice = 1;

shockmagnitude  = shockpars(1);
omega           = shockpars(3);

system_ext_assets = Tot_init_ext_assets;
system_deposits   = Tot_init_deposits;

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

    initshock_magnitude(i) = banks(i).shock(t);
end

tot_initshock_magnitude = sum(initshock_magnitude);


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
        % Updating balance interbank balance sheet entries for remaining active banks
        banks(i).assets.IB_Bil_Lending    = nonzeros(IB_loan_matrix(i,:,t));
        banks(i).assets.IB_Tot_Lending(t) = sum(banks(i).assets.IB_Bil_Lending);
        
        banks(i).liabilities.deposits(t+1)       = banks(i).liabilities.deposits(t); 
        banks(i).liabilities.IB_Bil_Borrowing    = nonzeros(IB_loan_matrix(:,i,t));
        banks(i).liabilities.IB_Tot_Borrowing(t) = sum(banks(i).liabilities.IB_Bil_Borrowing);
                         
        banks(i).status(t) = 1;
        banks(i).loss(t)   = 0;
        
        depositloss(i) = 0;
        
        banks(i).num_lenders(t+1)   = banks(i).num_lenders(t);
        banks(i).num_borrowers(t+1) = banks(i).num_borrowers(t);
            
    end
    
    for i = FailedBanks
        banks(i).status(t) = 0;
        banks(i).loss(t)   = 0;
    end

    shockedbanks_update = [];
    CurrentFailedBanks  = [];
    depositloss         = [];
    
    capitalloss_IBC     = []; % capital loss due to liquidity effects
    capitalloss_LE      = []; % capital loss due to interbank contagion
    
   if dispoutput == 1
        disp('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
        disp('-------------------------------- Shock absorption ------------------------------ ')
        disp('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
   end  
        
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
            
            banks(i).liabilities.capital(t+1)  = banks(i).liabilities.capital(t) - banks(i).shock(t); % Reduce capital by shock amount
            %banks(i).liabilities.deposits(t+1) = banks(i).liabilities.deposits(t);
            
            banks(i).loss(t) = banks(i).shock(t); 
            capitalloss_IBC(i) = banks(i).loss(t);
                        
            if dispoutput == 1
                fprintf(1,'Bank %d is able to absorb the shock in period %d with a capital loss of %.3f\n',i,t,banks(i).shock(t));
            end 
            
        % CASE II: Insufficient capital to absorb loss. Bank i fails and propagates to connected creditors 
        elseif banks(i).shock(t) >= banks(i).liabilities.capital(t) 
            
            failurecount = failurecount+1;
            
            CurrentFailedBanks(failurecount) = i;
 
            banks(i).status(t) = 0;
            banks(i).failtime  = t;
            banks(i).failed    = 1;
            
            % Firesale mechanism: Sell off all external assets
            banks(i).assets.externalassets(t+1) = 0;
            banks(i).firesales(t) = banks(i).assets.externalassets(t); 
                    
            shockedbanks_update = [shockedbanks_update banks(i).lender_ids']; % After i's failure, update set of shocked banks = i's creditors
            
            banks(i).liabilities.capital(t+1) = 0;
            
            banks(i).loss(t) = banks(i).liabilities.capital(t);
            capitalloss_IBC(i) = banks(i).loss(t);
            
            
            banks(i).residualshock(t) = banks(i).shock(t) - banks(i).liabilities.capital(t);
                        
            if dispoutput == 1
                if banks(i).num_lenders(t) ~= 0
                    fprintf(1,'Bank %d fails in period %d due to shock of size %.3f. Absorbs %.3f\n',...
                        i,t,banks(i).shock(t),banks(i).liabilities.capital(t));   
                else
                     fprintf(1,'Bank %d fails in period %d due to shock of size %.3f. No creditors for onward propagation\n',...
                        i,t,banks(i).shock(t));
                end
            end
            
            % Case II-1: Shock absorbed by interbank liabilities               
            if banks(i).residualshock(t)  < banks(i).liabilities.IB_Tot_Borrowing(t)
               
               banks(i).shocktransmit = ones(banks(i).num_lenders(t+1),1)*(1/banks(i).num_lenders(t+1)).*banks(i).residualshock(t);
               
               banks(i).liabilities.deposits(t+1) = banks(i).liabilities.deposits(t); % Customer deposits unaffected
                               
               if dispoutput == 1
                    fprintf(1,'--> Residual shock (%.3f) absorbed by interbank liabilities (%.3f). Transmit %.3f\n',...
                        banks(i).residualshock(t),sum(IB_loan_matrix(banks(i).lender_ids,i,t)),sum(banks(i).shocktransmit));
                    fprintf(1,'--> Updated interbank liabilities = %.3f\n',...
                        sum(IB_loan_matrix(banks(i).lender_ids,i,t)) - sum(banks(i).shocktransmit));
               end
                                                                 
            % Case II-2: Residual shock larger than capital AND interbank liabilities --> Customer deposits as ultimate sink
            else         
                banks(i).shocktransmit             = IB_loan_matrix(banks(i).lender_ids,i,t);
                banks(i).residualshock2(t)         = banks(i).residualshock(t) - banks(i).liabilities.IB_Tot_Borrowing(t);
                    
                if dispoutput == 1
                    fprintf(1,'--> Residual shock (%.3f) larger than interbank liabilities (%.3f). Transmit %.3f\n',...
                        banks(i).residualshock(t),sum(IB_loan_matrix(banks(i).lender_ids,i,t)),sum(banks(i).shocktransmit));
                     fprintf(1,'--> Updated interbank liabilities = %.3f\n',...
                       sum(IB_loan_matrix(banks(i).lender_ids,i,t)) - sum(banks(i).shocktransmit));
               end
                
                if banks(i).residualshock2(t)  < banks(i).liabilities.deposits(t) % Remaining residual absorbed by customer deposits
                    
                    banks(i).liabilities.deposits(t+1) = banks(i).liabilities.deposits(t) - banks(i).residualshock2(t);
                    
                    if dispoutput == 1
                    fprintf(1,'--> Sufficient deposits (%.3f) to absorb remaining residual shock (%.3f). Updated deposits = %.3f\n',...
                        banks(i).liabilities.deposits(t),banks(i).residualshock2(t),banks(i).liabilities.deposits(t+1));
                    end
                    
                    depositloss(i) = banks(i).liabilities.deposits(t) - banks(i).liabilities.deposits(t+1);
                    
                else % Customer deposits completely wiped out by remaining residual
                    
                    banks(i).liabilities.deposits(t+1) = 0;
                    
                    if dispoutput == 1
                    fprintf(1,'--> Insufficient deposits (%.3f) to absorbed remaining residual shock (%.3f). Updated deposits = 0\n',...
                        banks(i).liabilities.deposits(t),banks(i).liabilities.deposits(t+1));
                    end
                    
                    depositloss(i) = banks(i).liabilities.deposits(t);
                    
                end         
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
    end % Close current shocked bank for loop
    
    totaldepositloss(t)  = sum(depositloss);
    system_deposits(t+1) = system_deposits(t) - totaldepositloss(t);
            
    for i = 1:numel(shockedbanks)
        FSvec(i) = banks(shockedbanks(i)).firesales(t);
    end
    
    totalfiresales(t) = sum(FSvec);
    
    clearvars FSvec
    
    % Inverse demand function based on total external asset firesales
    
    if liquidityrisk == 1
        if totalfiresales(t) ~=0
            assetprice(t+1) = assetprice(t)*exp(-omega*(totalfiresales(t)/system_ext_assets(t)));
        else
            assetprice(t+1) = assetprice(t);
        end
    else
        assetprice(t+1) = assetprice(t);
    end
    
    system_ext_assets(t+1) = assetprice(t+1)*(system_ext_assets(t) - totalfiresales(t));

    % Output liquidity effect dynamics
    if dispoutput == 1 && assetprice(t+1) ~= assetprice(t) 
        disp('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
        disp('--------------------------- External asset dynamics ---------------------------')
        disp('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
        fprintf(1,'%.3f worth of external assets (= %.3f pc of total) were sold in firesale in period %d\n',...
            totalfiresales(t),(totalfiresales(t)/system_ext_assets(t))*100,t);
        fprintf(1,'--> Asset price change: %.3f -> %.3f\n',assetprice(t),assetprice(t+1));
        disp('------------------------------------------------------------------')
    elseif dispoutput == 1 && assetprice(t+1) == assetprice(t) 
        disp('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
        disp('--------------------------- External asset dynamics ---------------------------- ')
        disp('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
        fprintf(1,'No firesales in period %d\n',t)
        fprintf(1,'--> Asset price does not change from its previous value of %.3f\n',assetprice(t))
        disp('------------------------------------------------------------------')
    end
    
    % Output shock propagation dynamics now accounting for liquidity effects
    if dispoutput == 1
       disp('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
        disp('------------------------- Shock propagation dynamics -------------------------- ')
        disp('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
    end
    
    if ~isempty(CurrentFailedBanks)
     for i = CurrentFailedBanks
         
        Lvec_P = sprintf(' %.3g ',banks(i).lender_ids);
        
        if liquidityrisk == 1 % COND: Liquidity risk enabled and shocked bank fails in current period
            
            % Update shock transmission due to interbank contagion to allow for liquidity effects
            banks(i).ST_temp = banks(i).shocktransmit;
            banks(i).shocktransmit = banks(i).shocktransmit+...
                ones(banks(i).num_lenders(t+1),1)*(1/banks(i).num_lenders(t+1)).*(assetprice(t) - assetprice(t+1))*banks(i).firesales(t);
            
            IB_loan_matrix(banks(i).lender_ids,i,t+1)    = IB_loan_matrix(banks(i).lender_ids,i,t) - banks(i).shocktransmit;            
            IB_shockprop_matrix(banks(i).lender_ids,i,t) = banks(i).shocktransmit;
            
            if dispoutput == 1
                
                fprintf('Failed bank %d propagates %.3f to creditors (%.3f from shock and %.3f due to liquidity effects)\n',...
                     i,sum(banks(i).shocktransmit),sum(banks(i).ST_temp),sum(banks(i).shocktransmit)-sum(banks(i).ST_temp)); 
                 fprintf(1,'--> Firesales %.3f worth of external assets at current market price %.3f\n',...
                        banks(i).firesales(t),assetprice(t+1));
                 if banks(i).num_lenders(t) ~= 0
                    fprintf(1,'--> Transmits %.3f to each of %d creditors: [%s]\n',...
                        sum(banks(i).shocktransmit)/banks(i).num_lenders(t+1),banks(i).num_lenders(t+1),Lvec_P)
                 else
                     disp('no creditors to pass shock onto');
                 end
            end
        elseif liquidityrisk == 0
            if dispoutput == 1
                fprintf('Failed bank %d propagates %.3f to creditors\n',...
                     i,sum(banks(i).shocktransmit)); 
                 
                 if banks(i).num_lenders(t) ~= 0
                    fprintf(1,'--> Transmits %.3f to each of %d creditors: [%s]\n',...
                        sum(banks(i).shocktransmit)/banks(i).num_lenders(t+1),banks(i).num_lenders(t+1),Lvec_P)
                 else
                     disp('--> No creditors to pass shock onto');
                 end
            end               
        end
     end
    else
        if dispoutput == 1
            disp('Shock completely absorbed by remaining bank capital. No further propagation')
        end
    end
    
    if dispoutput == 1 && system_deposits(t+1) ~= system_deposits(t)
        fprintf(1,'System deposits updated from %.3f to %.3f\n',system_deposits(t),system_deposits(t+1));
    elseif dispoutput == 1 && system_deposits(t+1) == system_deposits(t)
        disp('Customer deposits unaffected in current round')
    end
     % Marking to market: Updating balance sheet to account for asset price changes in last period
     
     Fail_MTM_count = 0;
     Fail_MTM       = [];
         
     for i = ActiveBanks
         
     shockedbanks_update = unique(shockedbanks_update); % Control for the fact that failed banks can have creditors in common

        % Assets
        banks(i).assets.externalassets(t+1) = assetprice(t+1)*banks(i).assets.externalassets(t);
        if ~ismember(i,shockedbanks_update) % Active Banks not among current failed banks' creditors: 
            banks(i).assets.total(t+1)          = banks(i).assets.externalassets(t+1) + sum(IB_loan_matrix(i,:,t+1));
        else
            banks(i).assets.total(t+1)          = banks(i).assets.externalassets(t+1) + sum(IB_loan_matrix(i,:,t));
        end
            
        % Liabilities
        %banks(i).liabilities.deposits(t+1) = banks(i).liabilities.deposits(t);
        banks(i).liabilities.capital(t+1)  = banks(i).assets.total(t+1) - (sum(IB_loan_matrix(:,i,t+1)) +...
                                            banks(i).liabilities.deposits(t+1));                                                                                      
        banks(i).liabilities.total(t)      =   banks(i).assets.total(t); 
        
        %ea_init(i)   = banks(i).assets.externalassets(t);
        %ea_update(i) = banks(i).assets.externalassets(t+1);
     
        % Remove insolvent banks from the simulation
        if banks(i).liabilities.capital(t+1) < 0
            
            banks(i).failtime  = t;
            banks(i).failed    = 1;  
            banks(i).status(t) = 0; 
            banks(i).loss(t)   = banks(i).loss(t) + banks(i).liabilities.capital(t);
            
            capitalloss_LE(i)  = banks(i).liabilities.capital(t);
            
            Fail_MTM_count = Fail_MTM_count + 1;
            Fail_MTM(Fail_MTM_count) = i; % Vector containing banks who fail exclusively due to liquidity effects
            
            FailedBanks = [FailedBanks i];  
            ActiveBanks(ActiveBanks == i) = [];     
        end
     end
     
     % 2nd update of bank counterparty information
     for i = ActiveBanks
        for j = Fail_MTM
            
            banks(i).lender_ids(banks(i).lender_ids == j)   = [];
            banks(i).num_lenders(t+1) = numel(banks(i).lender_ids);
            
            banks(i).borrower_ids(banks(i).borrower_ids == j)   = [];
            banks(i).num_borrowers(t+1) = numel(banks(i).borrower_ids);
           
        end 
     end
          
     if dispoutput == 1
            disp('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
            disp('--------------------------- Mark to market accounting -------------------------- ')
            disp('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
            if assetprice(t+1) ~= assetprice(t)
                fprintf(1,'System external assets decrease from %.3f to %.3f due to liquidity effects\n',...
                     system_ext_assets(t), system_ext_assets(t+1));
                for i = Fail_MTM
                    fprintf(1,'Bank %d external assets revalued from %.3f to %.3f\n',...
                        i,banks(i).assets.externalassets(t),banks(i).assets.externalassets(t+1));
                    if ismember(i,shockedbanks_update)
                        disp('--> Creditor to a failed bank. Shock in next round')
                    else
                        disp('--> Not connected to any of the current failed banks')
                    end
                    fprintf(1,'--> Interbank assets updated from %.3f to %.3f\n',...
                        sum(IB_loan_matrix(i,:,t)),sum(IB_loan_matrix(i,:,t+1)));
                    fprintf(1,'--> Capital updated from %.3f to %.3f: INSOLVENT\n',...
                        banks(i).liabilities.capital(t),banks(i).liabilities.capital(t+1));
                end
                
                for i = setdiff(ActiveBanks,Fail_MTM)
                    fprintf(1,'Bank %d external assets revalued from %.3f to %.3f\n',...
                        i,banks(i).assets.externalassets(t),banks(i).assets.externalassets(t+1));
                    fprintf(1,'--> Capital updated from %.3f to %.3f: SOLVENT\n',...
                        banks(i).liabilities.capital(t),banks(i).liabilities.capital(t+1));             
                end
                
            elseif assetprice(t+1) == assetprice(t)
                fprintf(1,'No firesales in period %d therefore no change in value of system external assets (%.3f)\n',t, system_ext_assets(t));
                for i = Fail_MTM
                    fprintf(1,'Bank %d external assets revalued from %.3f to %.3f\n',...
                        i,banks(i).assets.externalassets(t),banks(i).assets.externalassets(t+1));
                    fprintf(1,'--> Interbank assets updated from %.3f to %.3f\n',...
                        sum(IB_loan_matrix(i,:,t)),sum(IB_loan_matrix(i,:,t+1)));
                 
                    fprintf(1,'--> Capital updated from %.3f to %.3f: INSOLVENT\n',...
                        banks(i).liabilities.capital(t),banks(i).liabilities.capital(t+1));
                end
            end
     end    
         
    if dispoutput == 1
        disp('Removal of insolvent banks from updated set of shocked banks')
        for j = Fail_MTM
            if ismember(j,shockedbanks_update)
            fprintf(1,'--> Bank %d declared insolvent due to asset price devaluation removed from updated set of shocked banks\n',j)
            end
        end
    end
    
    for j = FailedBanks
        shockedbanks_update(shockedbanks_update == j) = [];
    end
    
    NumFailedBanks(t) = failurecount+Fail_MTM_count; % Count total number of failed banks in current iteration

    totcapitalloss_IBC(t) = sum(capitalloss_IBC);
    totcapitalloss_LE(t)  = sum(capitalloss_LE);
    totcapitalloss(t)     = totcapitalloss_IBC(t) + totcapitalloss_LE(t);
    
    if dispoutput == 1
            disp('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
            disp('--------------------------------- Bank failures -------------------------------- ')
            disp('++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');
            fprintf(1,'%d banks failed in period %d\n',NumFailedBanks(t),t)
            fprintf(1,'--> %d due to interbank contagion and %d due to liquidity effects\n',failurecount,Fail_MTM_count);
            fprintf(1,'--> Total number of failed banks = %d\n',sum(NumFailedBanks));
            fprintf(1,'Total capital loss of %.3f in period %d\n',totcapitalloss(t),t)
            fprintf(1,'--> %.3f due to interbank contagion and %.3f due to liquidity effects\n',totcapitalloss_IBC(t),totcapitalloss_LE(t))
            fprintf(1,'--> Total capital loss = %.3f\n',sum(totcapitalloss))
    end
    
    if absorbcount == numel(shockedbanks) % Break loop condition 1: If shock is absorbed by all banks in current propagation round
        if dispoutput == 1
            disp('                       o----------------------------------o')
            fprintf('                       | Shock fully absorbed by period %d |\n',t);
            disp('                       o----------------------------------o')
        end
        shockterm  = t;
        assetprice = assetprice(1:shockterm+1); 
        break
    end
    
    if isempty(shockedbanks_update) % Break loop  condition 2: Current set of shocked banks do not have creditors for onward propagation
        if dispoutput == 1
        disp('                       o-------------------------------------------------------o')
            fprintf('                       | No more creditors for further propagation in period %d |\n',t);
        disp('                       o-------------------------------------------------------o')       
        end
        shockterm  = t;
        assetprice = assetprice(1:shockterm+1); 
        break
    end
    
    NumFailedBanks(t) = failurecount+Fail_MTM_count; % Count total number of failed banks in current iteration
    
    nonshockedbanks   = setdiff(ActiveBanks,shockedbanks);
    
    % Update values for banks that weren't shocked in the current period
    for i = nonshockedbanks
        %banks(i).liabilities.capital(t+1)   = banks(i).liabilities.capital(t);
        %banks(i).liabilities.deposits(t+1)  = banks(i).liabilities.deposits(t); 
        %banks(i).assets.externalassets(t+1) = banks(i).assets.externalassets(t);
        
        banks(i).status(t) = 1;
        banks(i).loss(t)   = 0;
    end
    shockedbanks = shockedbanks_update; % Update set of shocked banks to pass onto next iteration round
    
    t=t+1;
end % end while loop

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

TotDepositLoss      = sum(totaldepositloss);
Norm_TotDepositLoss = TotDepositLoss/Tot_init_deposits;

end