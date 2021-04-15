using Pkg
Pkg.add("StatsBase")
using StatsBase

function randomized_classical_shadow(num_total_measurements, system_size)
    # Implementation of the randomized classical shadow
    #
    #    num_total_measurements: int for the total number of measurement rounds
    #    system_size: int for how many qubits in the quantum system
    #
    measurement_procedure = [] #Initializes which measurements are used
    for t in 1:num_total_measurements
        single_round_measurement = [StatsBase.sample(["X", "Y", "Z"]) for i=1:system_size]
        push!(measurement_procedure, single_round_measurement)
    end
end

function derandomized_classical_shadow(all_observables, num_of_measurements_per_observable, system_size, weight=nothing)
   #
    # Implementation of the derandomized classical shadow
    #
    #     all_observables: a list of Pauli observables, each Pauli observable is a list of tuple
    #                      of the form ("X", position) or ("Y", position) or ("Z", position)
    #     num_of_measurements_per_observable: int for the number of measurement for each observable
    #     system_size: int for how many qubits in the quantum system
    #     weight: None or a list of coefficients for each observable
    #             None -- neglect this parameter
    #             a list -- modify the number of measurements for each observable by the corresponding weight
    #
    #if weight == nothing
    if weight == nothing 
        weight = ones(length(all_observables))
    end

    @assert length(weight)==length(all_observables)

    function cost_function(num_of_measurements_so_far, num_of_matches_needed_in_this_round)
        #Modified cost function in Huang2021 eq. C7/C8
        # num_of_measurements_so_far == number of measurements of "measurements per observable" covered in completely deterministic measurements
        # num_of_matches_needed_in_this_round == Per observable, how much would it need to match in the remaining qubits in this measurement i.e. weight of observable for qubit k to n

        eta = 0.9 # a hyperparameter subject to change
        nu = 1 - exp(-eta/2)

        cost = 0
        for (i, zipitem) in enumerate(zip(num_of_measurements_so_far, num_of_matches_needed_in_this_round))
            measurement_so_far, matches_needed = zipitem
            if num_of_measurements_so_far[i] >= floor(weight[i] * num_of_measurements_per_observable)
                continue #? Hard limit: wouldn't we want to maximize the cost for sufficiently measured observables?
            end

            if system_size < matches_needed
                V = eta / 2 * measurement_so_far
            else #Second term non-zero if all 1 to k qubit elements compatible 
                V = eta / 2 * measurement_so_far - log(1 - nu / (3^matches_needed))
            end
            cost += exp(-V/weight[i]) #mind the negative!
        end
        return cost
    end

    function match_up(qubit_i, dice_roll_pauli, single_observable)
        #Checks whether single qubit measurement at qubit i compatible with a single observable at qubit i
        #qubit_i: measurement qubit: int
        #dice_roll_pauli: Pauli guess: String
        #single_observable: single (full qubit) obsersable to check compatibility: tuple(str, int)


        for (pauli, pos) in single_observable #Searches through qubits for qubit i
            if pos != qubit_i
                continue
            else
                if pauli != dice_roll_pauli  #Compatible? NOTE: identities don't appear in single_observable *explicitely*!
                    return -1 #Incompatible (also completely)
                else
                    return 1 #Compatible on this qubit
                end
            end 
        end
        return 0
    end

    num_of_measurements_so_far = zeros(length(all_observables))
    measurement_procedure = []
    
    #No need to fix measurement number M: use modified cost function (M independent) -> also NO need to first randomize everything = Simulate dice rolls! See below.

    for repetition in 1:(num_of_measurements_per_observable * length(all_observables))
        #Iteration over measurements (if fixed = 1 to M)
        # A single round of parallel measurement over "system_size" number of qubits
        num_of_matches_needed_in_this_round = [length(P) for P in all_observables] # Recall, we have list of tuples = list of observables i.e. len(P) = observable weight
        single_round_measurement = []

        for qubit_i in 0:system_size-1 #Iteration over qubits within measurement m=repetition and gain cost
            cost_of_outcomes = Dict{String,Float64}("X" => 0, "Y" => 0, "Z" => 0)

            for dice_roll_pauli in ["X", "Y", "Z"] #Evaluate cost per Pauli option
                # Assume the dice rollout to be "dice_roll_pauli"
                
                for (i, single_observable) in enumerate(all_observables)
                    #Adjust number of matches needed: 
                    #if Pauli incompatible for a single observable -> let it blow up (so big that it also blows up for any proceeding iteration on the same measurement P )
                    #if compatible -> need to take reduction of LEFT accountable weight into sum of cost function of 2nd term in cost function
                    #
                    
                    result = match_up(qubit_i, dice_roll_pauli, single_observable)
                    if result == -1
                        num_of_matches_needed_in_this_round[i] += 100 * (system_size+10) # impossible to measure
                        #? Why this number: 100*(system_size+10)??, yes needs to be big enough, but why such specific? why not system_size+1?
                    end
                    if result == 1
                        num_of_matches_needed_in_this_round[i] -= 1 # m/atch up one Pauli X/YZ
                    end
                end
                
                cost_of_outcomes[dice_roll_pauli] = cost_function(num_of_measurements_so_far, num_of_matches_needed_in_this_round)

                # Revert the dice roll for next dice roll trail from same config
                for (i, single_observable) in enumerate(all_observables)
                    result = match_up(qubit_i, dice_roll_pauli, single_observable)
                    if result == -1
                        num_of_matches_needed_in_this_round[i] -= 100 * (system_size+10) # impossible to measure
                    end
                    if result == 1
                        num_of_matches_needed_in_this_round[i] += 1 # match up one Pauli X/Y/Z
                    end
                end
            end

            for dice_roll_pauli in ["X", "Y", "Z"] #Award minimal cost Pauli, store it and bookkeep number of measurements left within this measurement iteration m for all compatible observables 
                if minimum(collect(values(cost_of_outcomes))) < cost_of_outcomes[dice_roll_pauli] #Check next dice role if cost is not reduced by this dice role iteration
                    #Mind: if-expression checks if it isn't the minimal cost!
                    continue
                end
                
                # The best dice roll outcome will come to this line
                push!(single_round_measurement, dice_roll_pauli) #store
                for (i, single_observable) in enumerate(all_observables)
                    result = match_up(qubit_i, dice_roll_pauli, single_observable)
                    if result == -1
                        num_of_matches_needed_in_this_round[i] += 100 * (system_size+10) # impossible to measure
                    end
                    if result == 1
                        num_of_matches_needed_in_this_round[i] -= 1 # match up one Pauli X/Y/Z
                    end
                end
                break #First best Pauli that has smaller cost is chosen
            end
        end
        push!(measurement_procedure, single_round_measurement)

        for (i, single_observable) in enumerate(all_observables)
            if num_of_matches_needed_in_this_round[i] == 0 # finished measuring all qubits + observable i covered completely in this measurement?
                num_of_measurements_so_far[i] += 1
            end
        end

        success = 0
        for (i, single_observable) in enumerate(all_observables) #Check for how many observables the targeted goal of number of measurements has been fullfilled yet
            if num_of_measurements_so_far[i] >= floor(weight[i] * num_of_measurements_per_observable)
                success += 1
            end
        end

        status = convert(Int64,floor(num_of_measurements_per_observable*success/length(all_observables)))
        println("Measurement $repetition: reached $status measurements for all observables ")

        if success == length(all_observables)  #If we already reached our goal don't generate any new measurements i.e. enlarge M further
            break
        end
    end

    return measurement_procedure  #End and return measurement procedure, even if measurement number goal is not reached for all observables
end

#
# The following code is only used when we run this code through the command line interface
#

if abspath(PROGRAM_FILE) == @__FILE__
    function print_usage()
        println(stderr, "Usage:\n")
        println(stderr,"./shadow_data_acquisition -d [number of measurements per observable] [observable.txt]")
        println(stderr,"    This is the derandomized version of classical shadow.")
        println(stderr,"    We would output a list of Pauli measurements to measure all observables")
        println(stderr,"    in [observable.txt] for at least [number of measurements per observable] times.")
        println(stderr,"<or>\n")
        println(stderr,"./shadow_data_acquisition -r [number of total measurements] [system size]")
        println(stderr,"    This is the randomized version of classical shadow.")
        println(stderr,"    We would output a list of Pauli measurements for the given [system size]")
        println(stderr,"    with a total of [number of total measurements] repetitions.")
    end
    
    if length(ARGS) != 3
        print_usage()
    end

    if ARGS[1] == "-d"
        #open(ARGS[3], "r") do f
        #    content = f.readlines()
        #end
        content = readlines(ARGS[3])
        system_size = parse(Int64, content[1])

        all_observables = []
        weightArray = []

        for line in content[2:end]
            one_observable = []
            wBool = 0 #Does the provided file include weights?
            if (length(split(line, " "))-1) % parse(Int64, split(line, " ")[1]) != 0
                wBool = 1
            end

            for (pauli_XYZ, position) in zip(split(line, " ")[2:2:end-wBool], split(line," ")[3:2:end-wBool])
                push!(one_observable, (pauli_XYZ, parse(Int64, position)))
            end
            push!(all_observables, one_observable)

            if wBool == 1
                push!(weightArray, parse(Float64, split(line, " ")[end]))
            end
            
        end
        meas_repetitions_per_observable = parse(Int64, ARGS[2])
        if length(weightArray) == 0
            measurement_procedure = derandomized_classical_shadow(all_observables, meas_repetitions_per_observable, system_size)
        elseif length(weightArray) == length(all_observables)
            println("We detected weight inputs!")
            measurement_procedure = derandomized_classical_shadow(all_observables, meas_repetitions_per_observable, system_size, weightArray)
        else
            println(stderr,"Check your weights!")
            print_usage()
            measurement_procedure = []
        end
    elseif ARGS[1] == "-r"
        measurement_procedure = randomized_classical_shadow(parse(Int64, ARGS[2]), parse(Int64, ARGS[3]))
    else
        print_usage()
    end

    for single_round_measurement in measurement_procedure
        println(join(single_round_measurement, " "))
    end
end