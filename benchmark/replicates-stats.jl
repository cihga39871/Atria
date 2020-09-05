
using Statistics

if isempty(ARGS)
    println("""
    Usage: $(@__FILE__) FILE1 FILE2...

    Format of FILEs has to be the same. The numeric values in the same position of FILEs are replaced with their mean and standard deviation.
    """)
    exit()
end

contents = map(x -> read(x, String), ARGS)
cells = map(x -> split(x, r"[^\w\d\.\-\:]+"), contents)

cell_matrix = String.(hcat(cells...))

function Base.isnumeric(x::String)
    isfloat = !(tryparse(Float64, x) === nothing)
    istime = occursin(r"^(\d+:)+\d+(\.\d+)?$", x)
    isfloat | istime
end

function Base.isnumeric(x::AbstractArray)
    all(isnumeric.(x))
end

function parse_numeric(x::String)
    float = tryparse(Float64, x)
    if !isnothing(float)
        return float
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


mean_std_strings = map(eachrow(cell_matrix)) do vec
    if isnumeric(vec)
        vals = parse_numeric.(vec)
        std = Statistics.std(vals)
        if std != 0
            digit = -(floor(Int, log10(std)) - 1)
        else
            digit = 0
        end
        if digit <= 0
            std_string = round(Int, std) |> string
            mean_string = round(Int, Statistics.mean(vals)) |> string
        else
            std_string = round(std, sigdigits=2) |> string
            mean_string = round(Statistics.mean(vals), digits=digit) |> string
        end
        mean_string * " ± " * std_string
    else
        vec[1]
    end
end

specials = split(contents[1], r"[\w\d\.\-\:]+")

N = min(length(specials), length(mean_std_strings))
result = ""

if mean_std_strings[1] == ""
    for i = 1:N
        global result
        result *= mean_std_strings[i] * specials[i]
    end
else
    for i = 1:N
        global result
        result *= specials[i] * mean_std_strings[i]
    end
end

filename = "stats." * basename(ARGS[1])
write(filename, result)

@info "Done" Output=filename
