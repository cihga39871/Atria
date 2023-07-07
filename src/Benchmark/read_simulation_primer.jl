#!/usr/bin/env julia

# using ArgParse

using BioSequences

function parsing_args_simulate(args; exit_after_help = true)
    settings = ArgParseSettings(exit_after_help = exit_after_help)

    add_arg_group!(settings, "output")
    @add_arg_table! settings begin
        "--prefix", "-o"
            help = "prefix of output fastq files"
            metavar = "PREF"
            default = "read_simulation"
    end

    add_arg_group!(settings, "simulation")
    @add_arg_table! settings begin
        "--repeat", "-x"
            help = "repeat times for each case"
            default = 30000
            arg_type = Int64
        "--primer1", "-a"
            help = "read 1 primer"
            metavar = "SEQ"
            default = "AHCGATGAAGAACRYAG"
        "--primer2", "-A"
            help = "read 2 primer"
            metavar = "SEQ"
            default = "CTTATTGATATGCTTAAGTTCAG"
        "--seq-length", "-s"
            help = "a given sequence length; simulated sequence length might be 1 base more than the value because of simulated phasing error"
            default = 100
            arg_type = Int64
        "--insert-size-range", "-i"
            help = "range of insert size"
            nargs = '+'
            arg_type = Int64
            default = [80:2:120;]
        "--subsitution-rate", "-S"
            help = "subsitution rate per base. it is random for each base. error type includs mismatch"
            nargs = '+'
            arg_type = Float64
            default = [0.001:0.001:0.005;]
        "--insertion-rate", "-I"
            help = "insertion rate; number of arg should be the same as --subsitution-rate"
            nargs = '+'
            arg_type = Float64
            default = [1.0e-5:1.0e-5:5.0e-5;]
        "--deletion-rate", "-D"
            help = "deletion rate; number of arg should be the same as --subsitution-rate"
            nargs = '+'
            arg_type = Float64
            default = [1.0e-5:1.0e-5:5.0e-5;]
    end
    return parse_args(args, settings)
end


@inline function simulate_insert(insert_size::Int64)
    randdnaseq(insert_size)
end

"""
    simulate_error(base::DNA, sub_rate::Float64, insert_rate::Float64, del_rate::Float64)

Return `(base::DNA, iserror::Int64)`
"""
@inline function simulate_error(base::DNA, sub_rate::Float64, insert_rate::Float64, del_rate::Float64)
    bases = [DNA_A, DNA_T, DNA_C, DNA_G]

    randfloat = rand()
    if randfloat <= sub_rate
        ## subsitution
        idx_base = findfirst(x -> x == base, bases)
        idx_sub = rand(1:3)
        if idx_base == idx_sub
            return DNA_G, 1  # bases[4] == "G"
        else
            return bases[idx_sub], 1
        end
    else
        randfloat -= sub_rate
        if randfloat <= insert_rate
            ## insert
            res = base * rand(bases)
            return res, 1
        else
            randfloat -= insert_rate
            if randfloat <= del_rate
                ## deletion
                return DNA_Gap, 1
            end
            ## no error
            return base, 0
        end
    end
end

"""
    simulate_read(insert::String, primer::String, sub_rate::Float64, insert_rate::Float64, del_rate::Float64, seq_length::Int64)

Return `res, true_insert_size, nerror_insert, nerror_primer`
"""
@inline function simulate_read(insert::String, primer_head::LongDNA{4}, primer_tail::LongDNA{4}, sub_rate::Float64, insert_rate::Float64, del_rate::Float64, seq_length::Int64)
    res = LongDNA{4}
    ninsert = length(insert)
    nerror_insert = 0
    nerror_primer_head = 0
    nerror_primer_tail = 0





    # simulate insert
    the_base = ""
    for i in 1:ninsert
        the_base, iserror = simulate_error(insert[i:i], sub_rate::Float64, insert_rate::Float64, del_rate::Float64)
        n_the_base = length(the_base)
        res *= the_base
        nerror_insert += iserror
        current_insert_size = length(res)
        if current_insert_size >= seq_length
            if n_the_base <= 1
                return res, current_insert_size, nerror_insert, nerror_primer
            elseif n_the_base == 2  # error with insert. Inserted part is not of real DNA fragment!
                return res[1:end-1], current_insert_size - 1, nerror_insert, nerror_primer
            else
                @error "Bugs at simulate_read()" the_base iserror length(res) _module=nothing _group=nothing _id=nothing _file=nothing
            end
        end
    end
    if length(the_base) == 2  # error with insert. Inserted part is not of real DNA fragment!
        true_insert_size = length(res) - 1
    else
        true_insert_size = length(res)
    end

    nprimer = length(primer)
    for i in 1:nprimer
        the_base, iserror = simulate_error(primer[i:i], sub_rate::Float64, insert_rate::Float64, del_rate::Float64)
        res *= the_base
        nerror_primer += iserror
        if length(res) >= seq_length
            return res, true_insert_size, nerror_insert, nerror_primer
        end
    end

    nrandom = seq_length - length(res)
    res *= simulate_insert(nrandom)

    return res, true_insert_size, nerror_insert, nerror_primer
end

@inline function complement_char(c::Char)
    if c == 'A'
        'T'
    elseif c == 'T'
        'A'
    elseif c == 'C'
        'G'
    elseif c == 'G'
        'C'
    else
        'N'
    end
end

@inline function reverse_complement(s::String)
    ns = lastindex(s::String)
    char_vec = map(x -> complement_char(s[x]), ns:-1:1)
    string(char_vec...)
end

function writeseq(io::IO, header::String, seq::String; error_rate=0.0001)
    println(io, header)
    println(io, seq)
    println(io, "+")
    qual_char = if error_rate < 0.0001
        'J'
    else
        Char(round(Int, -10 * log10(error_rate)) + 33)
    end
    println(io, qual_char ^ length(seq))
end

function julia_wrapper_simulate(ARGS; exit_after_help = true)
    time0 = time()

    if length(ARGS) == 0
        parsing_args_simulate(["-h"], exit_after_help = exit_after_help)
        return 0
    end
    args = parsing_args_simulate(ARGS, exit_after_help = exit_after_help)
    args === nothing && return 0

    r1 = args["prefix"] * ".R1.fastq"
    r2 = args["prefix"] * ".R2.fastq"

    r1_io = open(r1, "w+")
    r2_io = open(r2, "w+")


    @info "read simulation: output files" r1 r2

    primer1 = args["primer1"]
    primer2 = args["primer2"]
    repeat_times = args["repeat"]
    seq_length = args["seq-length"]
    insert_sizes = args["insert-size-range"]

    insert_rates = args["insertion-rate"]
    deletion_rates = args["deletion-rate"]
    subsitution_rates = args["subsitution-rate"]

    length(insert_rates) == length(deletion_rates) == length(subsitution_rates) ||
        error("ArgumentError: the numbers of args of --subsitution-rate, --insertion-rate, and --deletion-rate should be the same. Abort.")

    error_rates = insert_rates .+ deletion_rates .+ subsitution_rates

    any(error_rates .> 1) &&
    error("ArgumentError: any dot sums of --subsitution-rate, --insertion-rate, and --deletion-rate should be less than one. Abort.")


    read_pair_count = repeat_times * length(insert_sizes) * length(error_rates)

    read_id = 0
    for insert_size in insert_sizes
        if insert_size >= 0
            for (i_rate, error_rate) in enumerate(error_rates)
                insert_rate = insert_rates[i_rate]
                del_rate = deletion_rates[i_rate]
                sub_rate = subsitution_rates[i_rate]

                for rep in 1:repeat_times
                    read_id += 1
                    insert = simulate_insert(insert_size::Int64)
                    r1_seq, r1_true_insert_size, r1_nerror_insert, r1_nerror_primer = simulate_read(insert, primer1, primer2_rc, sub_rate::Float64, insert_rate::Float64, del_rate::Float64, seq_length::Int64)
                    r2_seq, r2_true_insert_size, r2_nerror_insert, r2_nerror_primer = simulate_read(reverse_complement(insert), primer2, primer1_rc, sub_rate::Float64, insert_rate::Float64, del_rate::Float64, seq_length::Int64)

                    r1_header = "@PeReadSimulator2:$read_id:$rep TRUE=$r1_true_insert_size INSERT_SIZE=$insert_size ERROR_RATE=$error_rate SEQ_LENGTH=$seq_length ERROR_INSERT=$r1_nerror_insert ERROR_primer=$r1_nerror_primer SUB=$sub_rate INS=$insert_rate DEL=$del_rate"
                    r2_header = "@PeReadSimulator2:$read_id:$rep TRUE=$r2_true_insert_size INSERT_SIZE=$insert_size ERROR_RATE=$error_rate SEQ_LENGTH=$seq_length ERROR_INSERT=$r2_nerror_insert ERROR_primer=$r2_nerror_primer SUB=$sub_rate INS=$insert_rate DEL=$del_rate"

                    writeseq(r1_io, r1_header, r1_seq, error_rate=error_rate)
                    writeseq(r2_io, r2_header, r2_seq, error_rate=error_rate)
                end
            end
        end
    end
    close(r1_io)
    close(r2_io)
    @info "read simulation: all done" elapsed=time() - time0
    return 0
end
