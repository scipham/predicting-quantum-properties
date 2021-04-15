using Pkg

function hammingDistance(n1, n2)
    x = n1 ⊻ n2 
    setBits = 0
  
    while x > 0
        setBits += x & 1
        x >>= 1
    end
    return setBits
end

function alt_predict_renyi_2_entropy(full_measurement, all_subsystems)

    for (s,subsystem) in enumerate(all_subsystems)
        subsystem_size = length(subsystem)

        measurements_covered = []
        measurements_covered_outcomes = []
        unzip(a) = [getindex.(a, i) for i in 1:length(a[1])]
    
        for (m,single_measurement) in enumerate(full_measurement)
            #cnt = count(x-> x==single_measurement, full_measurement)
            sm_operator, sm_outcome = unzip(single_measurement)
            push!(measurements_covered, sm_operator)
            push!(measurements_covered_outcomes, sm_outcome)
        end

        unique_measurements = unique(measurements_covered)
        unique_m_bitstring_prob = []
 

        for (m, single_unique_measurement) in enumerate(unique_measurements)
            @show m
            unique_idxs = findall(x-> x==single_unique_measurement, measurements_covered)
            
            N_M = length(unique_idxs)

            unique_m_bitstrings = []
            for i in unique_idxs
                ittobit = convert(Array{Int64},(measurements_covered_outcomes[i].+1)./2) #-/+1 -> 0,1
                push!(unique_m_bitstrings, [ittobit[subsystem[si]+1] for si in 1:subsystem_size])
            end
            
            bitstring_sum_m = zeros(2^subsystem_size)

            for umb in unique_m_bitstrings
                found = 0
                umb = string.(umb)
                pushfirst!(umb, "0b") #prepare binary-to-int conversion
                for b in 0:(2^(subsystem_size)-1)
                    if b == parse(Int64, join(umb))
                        bitstring_sum_m[b+1] += 1
                        found = 1
                    end
                end
                if found == 0
                    @show "Oops, something goes wrong"
                end
            end

            bitstring_prob_m = bitstring_sum_m ./ length(unique_m_bitstrings)
            push!(unique_m_bitstring_prob, bitstring_prob_m)
        end

        totalProb = 0
        sum_terms_ij = []
        for i in 0:(2^(subsystem_size) -1)
            for j in 0:(2^(subsystem_size) -1)
                P_ij_mn_sum = 0
                for (m, m_prob_m) in enumerate(unique_m_bitstring_prob)
                    for (n, m_prob_n) in enumerate(unique_m_bitstring_prob)
                        @show m,n
                        if m != n
                            P_ij_mn_sum += m_prob_m[i+1] * m_prob_n[j+1]
                        end
                    end
                end
                P_ij_mean = P_ij_mn_sum/(length(unique_m_bitstring_prob)*(length(unique_m_bitstring_prob)-1))
                totalProb += P_ij_mean
                push!(sum_terms_ij, (-2)^(-convert(Float64,hammingDistance(i,j))) * P_ij_mean)
                #@show i,j, P_ij_mn_sum/(length(unique_m_bitstring_prob)*(length(unique_m_bitstring_prob)-1))
            end
        end
        @show totalProb
        trace_of_square = 2^subsystem_size * sum(sum_terms_ij)
        println(-1.0 * log2(min(max(trace_of_square, 1.0 / (2.0^subsystem_size)), 1.0 - 1e-9)))
        @show "ending on subsystem"
    end
end

function alt2_predict_renyi_2_entropy(full_measurement, all_subsystems)

    for (s,subsystem) in enumerate(all_subsystems)
        subsystem_size = length(subsystem)

        measurements_covered = []
        measurements_covered_outcomes = []
        unzip(a) = [getindex.(a, i) for i in 1:length(a[1])]
    
        for (m,single_measurement) in enumerate(full_measurement)
            #cnt = count(x-> x==single_measurement, full_measurement)
            sm_operator, sm_outcome = unzip(single_measurement)
            push!(measurements_covered, [sm_operator[si+1] for si in 1:subsystem_size])
            push!(measurements_covered_outcomes, [sm_outcome[si+1] for si in 1:subsystem_size])
        end

        unique_measurements = unique(measurements_covered)
        unique_m_bitstring_prob = []
 

        for (m, single_unique_measurement) in enumerate(unique_measurements)
            unique_idxs = findall(x-> x==single_unique_measurement, measurements_covered)
            
            N_M = length(unique_idxs)

            unique_m_bitstrings = []
            for i in unique_idxs
                push!(unique_m_bitstrings, convert(Array{Int64},(measurements_covered_outcomes[i].+1)./2))
            end
            
            bitstring_sum_m = zeros(2^subsystem_size)

            for umb in unique_m_bitstrings
                found = 0
                umb = string.(umb)
                pushfirst!(umb, "0b") #prepare binary-to-int conversion
                for b in 0:(2^(subsystem_size)-1)
                    if b == parse(Int64, join(umb))
                        bitstring_sum_m[b+1] += 1
                        found = 1
                    end
                end
                if found == 0
                    @show "Oops, something goes wrong"
                end
            end

            bitstring_prob_m = bitstring_sum_m ./ length(unique_m_bitstrings)
            push!(unique_m_bitstring_prob, bitstring_prob_m)
        end

        totalProb = 0
        sum_terms_ij = []
        for i in 0:(2^(subsystem_size) -1)
            for j in 0:(2^(subsystem_size) -1)
                P_ij_mn_sum = 0
                for (m, m_prob_m) in enumerate(unique_m_bitstring_prob)
                    for (n, m_prob_n) in enumerate(unique_m_bitstring_prob)
                        if m != n
                            P_ij_mn_sum += m_prob_m[i+1] * m_prob_n[j+1]
                        end
                    end
                end
                P_ij_mean = P_ij_mn_sum/(length(unique_m_bitstring_prob)*(length(unique_m_bitstring_prob)-1))
                totalProb += P_ij_mean
                push!(sum_terms_ij, (-2)^(-convert(Float64,hammingDistance(i,j))) * P_ij_mean)
                #@show i,j, P_ij_mn_sum/(length(unique_m_bitstring_prob)*(length(unique_m_bitstring_prob)-1))
            end
        end
        @show totalProb
        trace_of_square = 2^subsystem_size * sum(sum_terms_ij)
        println(-1.0 * log2(min(max(trace_of_square, 1.0 / (2.0^subsystem_size)), 1.0 - 1e-9)))
        @show "ending on subsystem"
    end
end

#Predicting Entanglement entropy: Experimental!
mapPauli = Dict("X" => 0, "Y"=> 1, "Z"=> 2)

function predict_renyi_2_entropy(full_measurement, all_subsystems)
    renyi_sum_of_binary_outcome = []
    renyi_number_of_outcomes = []

    for (s,subsystem) in enumerate(all_subsystems)
        subsystem_size = length(subsystem)

        renyi_sum_of_binary_outcome = zeros(2^(2*subsystem_size))
        renyi_number_of_outcomes = zeros(2^(2*subsystem_size))

        for (m, measurement) in enumerate(full_measurement)
            encoding = 0+1
            cumulative_outcome = 1
            
            renyi_sum_of_binary_outcome[1] += 1
            renyi_number_of_outcomes[1] += 1

            #Using gray code iteration over all 2^n possible outcomes
            for b in 1:(2^subsystem_size - 1) #Iterate over binary bit strings
                change_i = trailing_zeros(b)+1  #Look for next gray code flip bit from binary string
                index_in_original_system = subsystem[change_i]+1  #In the larger system, look up what this corresponds to
                cumulative_outcome *= measurement[index_in_original_system][2]  #[2] = outcome, record system qubit outcome for total
                @show change_i
                @show index_in_original_system
                @show mapPauli[measurement[index_in_original_system][1]]
                @show encoding
                encoding ⊻= (mapPauli[measurement[index_in_original_system][1]]+1) << (2*change_i)  #[1] = pauli X/Y/Z
                #Shifts Pauli operator number bitwise according to change_i occurence -> XOR with current bitstring = only change a single bit if operator AND shift equal to before = same encoding
                @show encoding
                renyi_sum_of_binary_outcome[encoding] += cumulative_outcome #Record total expectation (sample) sum for operator on shift encoding
                renyi_number_of_outcomes[encoding] += 1 #For probability normalization, count encoding coincidences
            end
        end 

        level_cnt=zeros(2 * subsystem_size) #level_count
        level_ttl=zeros(2 * subsystem_size) #level_total

        for c in 0:(2^(2*subsystem_size) -1) #Iterate over binary bit-strings
            nonId = 1
            for i in 1:subsystem_size
                nonId += convert(Int64, ((c>>(2*i)) & 3) != 0)  #bool -> 0 or 1 
                #? Undo encoding for subsystem qubit (=change_i) and bitstring; check whether 3=Z is involved in this bitstring
                #-> if it is add -> nonId = count for how many subsystem qubits this holds for this bit string c
            end
            if renyi_number_of_outcomes[c+1] >= 2
                level_cnt[nonId] += 1 #Record for this bitstring if number of measurements is significant (>=2) and record for which number of qubit Z-pauli condition applied to
            end
            level_ttl[nonId] += 1 #Record how many bitstrings considered for this for (probability) normalization later
        end

        trace_of_square = 0 #Trace of square of subsystem density matrix in prediction
        for c in 0:(2^(2*subsystem_size) -1)
            if renyi_number_of_outcomes[c+1] <= 1
                continue  #Insignificant/insufficient amount of measurement recordings for this bitstring
            end
            nonId = 1
            for i in 1:subsystem_size
                nonId += convert(Int64, ((c>>(2*i)) & 3) != 0)  #bool -> 0 or 1 
                #see comment above
            end
           trace_of_square += 1.0/(renyi_number_of_outcomes[c+1] * (renyi_number_of_outcomes[c+1] - 1)) * (renyi_sum_of_binary_outcome[c+1] * renyi_sum_of_binary_outcome[c+1] - renyi_number_of_outcomes[c+1])/(2^subsystem_size) * level_ttl[nonId] / level_cnt[nonId]
           #@show renyi_sum_of_binary_outcome[c+1] 
           #@show renyi_number_of_outcomes[c+1] 
           #@show level_ttl[nonId] 
           #@show level_cnt[nonId]
        end

        println(-1.0 * log2(min(max(trace_of_square, 1.0 / (2.0^subsystem_size)), 1.0 - 1e-9)))
    end
end

#Need "observables working on qubit i"-array (indices)! observables_acting_on_ith_qubit
#Nested? i.e. position -> which of the three paulis -> individual count 
#observables_acting_on_ith_qubit




#This is the real stuff:
function estimate_exp(full_measurement, one_observable)
    #Takes full measurement scheme with a single observable to predict its expectation
    sum_product, cnt_match = 0, 0 #Initialize 
    for single_measurement in full_measurement
        not_match = 0
        product = 1

        for (pauli_XYZ, position) in one_observable #Iterate over each observable qubit-part 
            if pauli_XYZ != single_measurement[position+1][1] #Check if compatible (identities, again, not included!)
                not_match = 1 
                break
            end
            product *= single_measurement[position+1][2]  #If compatible, "add" result to total expectation product
        end

        if not_match == 1
            continue
        end

        sum_product += product  #Add single observable expectation for each comptible measurement together (later taken average with)
        cnt_match += 1  #Count how many measurements were compatible for taking expectation average (/normalize sum) over all compatible measurements
    end
    return sum_product, cnt_match
end

if abspath(PROGRAM_FILE) == @__FILE__
    function print_usage()
        print(stderr, "Usage:\n")
        print(stderr, "./prediction_shadow -o [measurement.txt] [observable.txt]")
        print(stderr, "    This option predicts the expectation of local observables.")
        print(stderr, "    We would output the predicted value for each local observable given in [observable.txt]")
    end
    @show ARGS
    
    if length(ARGS) != 3
        print_usage()
        print("\n")
        exit()
    end

    measurements = readlines(ARGS[2])
    system_size = parse(Int64, measurements[1])
    
    full_measurement = []
    for line in measurements[2:end]
        single_measurement = []
        for (pauli_XYZ, outcome) in zip(split(line, " ")[1:2:end], split(line, " ")[2:2:end])
            push!(single_measurement, (pauli_XYZ, parse(Int64, outcome)))
        end
        push!(full_measurement, single_measurement)
    end

    if ARGS[1]=="-o"
        observablesContent = readlines(ARGS[3])
        @assert system_size == parse(Int64, observablesContent[1])

        weightArray = []
        
        for (obs_i, line) in enumerate(observablesContent[2:end])
            one_observable = []

            wBool = 0 #Does the provided file include weights?
            if (length(split(line, " "))-1) % parse(Int64, split(line, " ")[1]) != 0
                wBool = 1
            end

            for (pauli_XYZ, position) in zip(split(line," ")[2:2:end-wBool], split(line," ")[3:2:end-wBool])
                push!(one_observable, (pauli_XYZ, parse(Int64, position)))
            end

            if wBool == 1
                push!(weightArray, parse(Float64, split(line, " ")[end]))
            end

            sum_product, cnt_match = estimate_exp(full_measurement, one_observable)
            if cnt_match > 0
                predictedExp = sum_product / cnt_match #Final expectation prediction: Sum of average outcome in each measurement -> normalize = average
            else
                println("WARNING: Observable $obs_i not measured once, we'll expect 0 for it!")
                predictedExp = 0
            end
            println("For observable $obs_i : $predictedExp")
        end
    elseif ARGS[1]=="-e"
        subsystemsContent = readlines(ARGS[3]) #Read lines of subsystems file
        
        @assert system_size == parse(Int64, subsystemsContent[1])
        
        all_subsystems = []
        subsystem_sizes = []

        for line in subsystemsContent[2:end]
            one_subsystem = []
            subsystem_size = parse(Int64, split(line, " ")[1])
            for (i, qubit_i) in enumerate(split(line," ")[2:end])
                push!(one_subsystem, parse(Int64, qubit_i))
            end
            push!(subsystem_sizes, subsystem_size)
            push!(all_subsystems, one_subsystem)
        end
        alt_predict_renyi_2_entropy(full_measurement, all_subsystems)
    else
        println(stderr, "Check you arguments! \n")
        print_usage()
        println(stderr, "\n")
        exit()
    end
end