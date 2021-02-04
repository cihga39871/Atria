#!julia

using DataFrames, CSV

if isempty(ARGS) || !isfile(ARGS[1])
    println("""
    Usage: $(@__FILE__) time_benchmark.txt num_bases [stderr.pigz.log]

    time_benchmark.txt is the result of `pasteTimeOutput` (see simulate-run-bench.bash);
    num_bases is the number of bases processed.
    stderr.pigz.log compensate the bug of GNU TIME which cannot stat the subprocess (pigz) of Julia. The file is the result of ```
        /usr/bin/time pigz -p 8 -c Atria-consensus/*atria.fq 1>/dev/null 2>> stderr.pigz.log
        /usr/bin/time pigz -p 8 -c Atria/*atria.fq 1>/dev/null 2>> stderr.pigz.log
        /usr/bin/time pigz -cd $r1 $r2 > /dev/null 2>> stderr.pigz.log
    ```

    Result output to stdout.
    """)
    exit()
end

function parse_numeric(x::String)
    float = tryparse(Float64, x)
    if !isnothing(float)
        return float
    end
    # parse percentage
    if occursin(r"\%$", x)
        float = tryparse(Float64, x[1:end-1])
        if !isnothing(float)
            return float/100
        end
    end
    # parse time as D:H:M:S.MS
    xs = split(x, ":") |> reverse!
    second = parse(Float64, xs[1])
    for i in 2:length(xs)
        if i == 2
            second += 60 * parse(Float64, xs[i])
        elseif i == 3
            second += 3600 * parse(Float64, xs[i])
        elseif i == 4
            second += 24 * 3600 * parse(Float64, xs[i])
        else
            error("Failed to parse $x as the time format D:H:M:S")
        end
    end
    second
end

const THREADS_STR = ["-threads", "-thread", "-cores", "-t"]
function get_threads(x; THREADS_STR=THREADS_STR)
    thread = 1
    for thread_str in THREADS_STR
        m = match(Regex("$thread_str[= ]*([\\d]+)"), x)
        isnothing(m) && continue
        if length(m.captures) == 1
            thread = parse(Int, m.captures[1])
        end
    end
    @warn "Set threads == $thread for command: $x"
    thread
end

df = CSV.File(ARGS[1], header=false) |> DataFrame!
NUM_BASES = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1

const USERTIME_STR = "User time (seconds): "
const SYSTEMTIME_STR = "System time (seconds): "
const CPU_STR = "Percent of CPU this job got: "
const ELAPSEDTIME_STR = "Elapsed (wall clock) time (h:mm:ss or m:ss): "
const MEMORY_STR = "Maximum resident set size (kbytes): "

USERTIME = findfirst(x -> typeof(x)<:AbstractString && occursin(USERTIME_STR, x), df[1,1:end])
SYSTEMTIME = findfirst(x -> typeof(x)<:AbstractString && occursin(SYSTEMTIME_STR, x), df[1,1:end])
CPU = findfirst(x -> typeof(x)<:AbstractString && occursin(CPU_STR, x), df[1,1:end])
ELAPSEDTIME = findfirst(x -> typeof(x)<:AbstractString && occursin(ELAPSEDTIME_STR, x), df[1,1:end])
MEMORY = findfirst(x -> typeof(x)<:AbstractString && occursin(MEMORY_STR, x), df[1,1:end])

usertimes = parse_numeric.(replace.(df[!, USERTIME], USERTIME_STR=>""))
systemtimes = parse_numeric.(replace.(df[!, SYSTEMTIME], SYSTEMTIME_STR=>""))
# cpus = parse_numeric.(replace.(df[!, CPU], CPU_STR=>""))
elapsedtimes = parse_numeric.(replace.(df[!, ELAPSEDTIME], ELAPSEDTIME_STR=>""))
memories = parse_numeric.(replace.(df[!, MEMORY], MEMORY_STR=>""))
threads = get_threads.(df[!,1])

if length(ARGS) == 3
    #=
    stderr.pigz.log compensate the bug of GNU TIME which cannot stat the subprocess (pigz) of Julia. The file is the result of ```
        /usr/bin/time pigz -p 8 -c Atria-consensus/*atria.fq 1>/dev/null 2>> stderr.pigz.log
        /usr/bin/time pigz -p 8 -c Atria/*atria.fq 1>/dev/null 2>> stderr.pigz.log
        /usr/bin/time pigz -cd $r1 $r2 > /dev/null 2>> stderr.pigz.log
    ```
    =#
    pigz_time_file = ARGS[3]
    usertimes_pigz = parse_numeric.(readlines(pipeline(`grep -oE "[0-9\.\:]+user" $pigz_time_file`, `sed 's/user//'`)))
    systemtimes_pigz = parse_numeric.(readlines(pipeline(`grep -oE "[0-9\.\:]+system" $pigz_time_file`, `sed 's/system//'`)))

    rows_atria = map(x -> occursin(r"atria", x), df[!,1])
    rows_atria_no_consensus = map(x -> occursin(r"atria .*--no-consensus", x), df[!,1])
    rows_atria_consensus = rows_atria .⊻ rows_atria_no_consensus

    # add decompressing time
    usertimes[rows_atria] .+= usertimes_pigz[3]
    systemtimes[rows_atria] .+= systemtimes_pigz[3]

    # add compressing time
    usertimes[rows_atria_no_consensus] .+= usertimes_pigz[2]
    systemtimes[rows_atria_no_consensus] .+= systemtimes_pigz[2]

    usertimes[rows_atria_consensus] .+= usertimes_pigz[1]
    systemtimes[rows_atria_consensus] .+= systemtimes_pigz[1]
end

cpus = @. (usertimes + systemtimes) / elapsedtimes
efficiencies = @. NUM_BASES / elapsedtimes / cpus / 10^6 # M Bases/s/CPU
speeds = NUM_BASES ./ elapsedtimes / 10^6 # M Bases/s

dfout = DataFrame(
    "Threads" => threads,
    "Command" => df[!,1],
    "Efficiency (M Bases/s/CPU)" => efficiencies,
    "Speed (M Bases/s)" => speeds,
    "UserTime (s)" => usertimes,
    "SystemTime (s)" => systemtimes,
    "CPU" => cpus,
    "ElapsedTime (s)" => elapsedtimes,
    "MaxMemory (kB)" => memories
)

sort!(dfout, :Threads)

CSV.write(stdout, dfout; delim='\t')
