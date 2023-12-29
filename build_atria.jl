#!julia --color=yes

using Pkg
using UUIDs

if !isempty(intersect(ARGS, ["-h", "--help", "--h"]))
    println("Usage: $(@__FILE__) [path-to-install]")
    exit(1)
end

function check_compatibility()
    if Sys.iswindows()
        # without pre-compilation
        error("Atria does not support Windows Platform.")
        exit()
    end

    if !(v"1.8" <= VERSION < v"1.10")
        @warn "Julia version is not v1.8 or 1.9. The build might be fail." JULIA_VERSION = VERSION
    end

    if !(v"1.8" <= VERSION < v"1.9")
        @warn "It is recommended to build Atria using Julia v1.8 because it is 3-20% faster than v1.9." JULIA_VERSION = VERSION
    end

    try
        run(pipeline(`pigz --version`, stderr=devnull))
    catch
        error("Dependency missing: pigz not found.")
        exit()
    end
    try
        run(pipeline(`pbzip2 --version`, stderr=devnull))
    catch
        error("Dependency missing: pbzip2 not found.")
        exit()
    end

    # Pkg.add("PackageCompiler")
end

function check_biosequences_for_atria()
    biosequences_info = Pkg.dependencies()[UUID("7e6ae17a-c86d-528c-b3b9-7f778a29fe59")]
    if !occursin("cihga39871/BioSequences.jl", biosequences_info.git_source)
        # Found package BioSequences from $(biosequences_info.git_source). Atria needs a specific version of BioSequences.
        Pkg.add(url="https://github.com/cihga39871/BioSequences.jl")
    end
end

check_compatibility()

cd(@__DIR__)
Pkg.activate(".")
Pkg.update()
Pkg.add(url="https://github.com/cihga39871/BioSequences.jl")
check_biosequences_for_atria()
Pkg.resolve()
Pkg.instantiate()
Pkg.precompile()



import Pkg; Pkg.add("PackageCompiler")
using PackageCompiler
using Dates

ver = "atria"
if isfile("Project.toml")
    for i in readlines("Project.toml")
        if occursin(r"^version", i)
            ind = findall("\"", i)
            global ver = i[ind[1].start + 1 : ind[2].start - 1]
            break
        end
    end
end

if length(ARGS) == 1
    app_path = joinpath(ARGS[1], "app-$ver")
else
    app_path = joinpath(".", "app-$ver")
end

if isdir(app_path)
    app_path *= "_" * replace(string(now()), r":\d\d\..*" => "", ":" => "-")
end

bash_wrapper = raw"""
#!/usr/bin/env bash

SOURCE="${BASH_SOURCE[0]}"
while [ -h "$SOURCE" ]; do # resolve $SOURCE until the file is no longer a symlink
  DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"
  SOURCE="$(readlink "$SOURCE")"
  [[ $SOURCE != /* ]] && SOURCE="$DIR/$SOURCE" # if $SOURCE was a relative symlink, we need to resolve it relative to the path where the symlink file was located
done
DIR="$( cd -P "$( dirname "$SOURCE" )" >/dev/null 2>&1 && pwd )"

ALLARGS=( "$@" )

# default
JULIA_NUM_THREADS=8

while [[ $# -gt 0 ]]
do
    case "$1" in
        --threads)
        JULIA_NUM_THREADS="$2"
        shift
        shift
        ;;
        -t)
        JULIA_NUM_THREADS="$2"
        shift
        shift
        ;;
        *)
        shift
        ;;
    esac
done

# "$DIR/AtriaEntry" "${ALLARGS[@]}"

unameOut="$(uname -s)"
case "${unameOut}" in
    Linux*)     EXT=so;;
    Darwin*)    EXT=dylib;;
    *)          EXT=so;;
esac

ATRIA="$DIR/_atria"

"$ATRIA" "${ALLARGS[@]}" --julia-args --check-bounds=no --math-mode=fast -t "$JULIA_NUM_THREADS"

"""

precompile_execution_file = joinpath("test", "runtests.jl")

create_app(".", app_path, incremental = true, force = true, filter_stdlibs = false, sysimage_build_args = `-O3 --check-bounds=no --math-mode=fast`, precompile_execution_file = precompile_execution_file, executables = ["_atria" => "julia_main"])

# ext = Sys.isapple() ? "dylib" : "so"
# isfile(joinpath(app_path, "bin", "AtriaEntry.$ext"))

# copy julia exe
# julia_exe = joinpath(Sys.BINDIR, Base.julia_exename())
# cp(julia_exe, joinpath(app_path, "bin", Base.julia_exename()))

Atria_bash = joinpath(app_path, "bin", "atria")
io = open(Atria_bash, "w+")
write(io, bash_wrapper)
close(io)
chmod(Atria_bash, 0o755)

@info "Success. Atria is installed at $app_path/bin/atria"

check_compatibility()
