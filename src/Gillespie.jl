 #module MultilevelGenesis

#############################################
###### Dependencies #########################
#############################################
using Random
#using Plots
using DataFrames
using JSON
using CSV
using Combinatorics


#############################################
###### Rection Functions ####################
#############################################

function CalculatePropensities(concentrations, CRS) 
    ########## NOT UNIT TESTED! ##########
    PropensityVec = zeros(length(CRS.reaction_list))
    for rxn in CRS.reaction_list
        Ap = 0.0
        if rxn.propensity == "STD"
            Ap = STD_propensity(rxn, concentrations)
        end
        PropensityVec[rxn.index] = Ap
    end
    
    return PropensityVec
end

function PickReactionID(PropensityVec)
    Ap_tot = sum(PropensityVec)
    dice_roll = Ap_tot*rand()
    checkpoint = 0.0
    
    rID = 0
    while checkpoint < dice_roll
        rID += 1
        checkpoint += PropensityVec[rID]         
    end
   return rID 
end

function ExecuteReaction(concentrations, CRS, rID)
    ########## NOT UNIT TESTED! ##########
    rxn = CRS.reaction_list[rID]
    for (i, r) in enumerate(rxn.reactants)
        concentrations[r] -= rxn.react_coef[i]
    end
    for (j, p) in enumerate(rxn.products)
        concentrations[p] += rxn.prod_coef[j]
    end
    if any(x->x<0, concentrations)
        error("Concentrations went negative on reaction $rID")
    end
    
    return concentrations
end

function STD_propensity(rxn, concentrations)
    ########## NOT UNIT TESTED! ##########
    Ap  = rxn.rate_constant
    for (ri,r) in enumerate(rxn.reactants)
        rconc = concentrations[r]
        r_order = rxn.react_coef[ri]
        Ap = Ap*(prod([(rconc-i) for i in 0:(r_order-1)])) 
    end
    cat_enhance = 0.0
    for (ci,c) in enumerate(rxn.catalysts)
        cconc = concentrations[c]
        cat_enhance += cconc*rxn.cat_coef[ci]
    end
    Ap = Ap*(1 + cat_enhance)
    return Ap
end



#############################################
##### Stochastic Evolution Functions ########
#############################################

 
function Run_Binary_Ligation(nA::Int64, nB::Int64, max_L::Int64, kf::Float64, kb::Float64, volume::Float64, tau_max::Float64, tau_freq::Float64 = 0.1, T::Float64 = 300.0, R::Float64= 8.3144598, mA::Int64 = 100, mB::Int64 = 100)
    ### Unique Experiment number required
    
    ## Check for the correct output files and directories
    param_file = "../data/Binary_Ligation_Run_Parameters.csv"
    save_dir = "../data/Binary_Ligation_timeseries"
    
    if !isdir(save_dir)
        mkdir(save_dir)
    end
    
        
    ## Initialize Time values
    tau = 0.0
    freq_count = 0.0
    
    binary_CRS = generate_binary_ligation_CRS(max_L, kf, kb, volume, T)
    ### Redo output initialization using DataFrames
    output_DF = InitializeOutput(binary_CRS)
               
    nmolecules = length(binary_CRS.molecule_list)
    concentrations= zeros(nmolecules)
    
    Aindex = binary_CRS.molecule_dict["A"]
    Bindex = binary_CRS.molecule_dict["B"]
    concentrations[Aindex] = nA
    concentrations[Bindex] = nB
        
    #Initialize Propensities
    propensities = CalculatePropensities(concentrations, binary_CRS)

    ####### Main LOOP #######
    while tau < tau_max
        # Pick Reaction
        rID = PickReactionID(propensities)
        # Execute Reaction
        concentrations = ExecuteReaction(concentrations, binary_CRS, rID)
        # Calculate Propensities
        propensities = CalculatePropensities(concentrations, binary_CRS)
        # Record Data
        if tau >= freq_count
            output_DF[Symbol(freq_count)] = concentrations[:] 
            freq_count += tau_freq
            println(tau)
        end

        #Update Time
        Ap_tot = sum(propensities)
        tau -= (log(rand())/Ap_tot) 
    end

    # Save time series and params
    save_name, parameter_df = generate_output_data(binary_CRS, save_dir)#
    println(parameter_df)
    if !isfile(param_file)
        CSV.write(param_file, parameter_df)
    else
        runs_df = CSV.read(param_file)
        runs_df = [runs_df; parameter_df]
    end
    
    CSV.write(save_name, output_DF)
    return concentrations
end
#############################################
#end #### module #############################
#############################################