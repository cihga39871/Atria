#!julia --color=yes

using Pkg

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

    if !(v"1.4" <= VERSION < v"1.5")
        @error "Performance: Julia version is not v1.4. The build might be fail. Atria built with Julia v1.5.1 is slower than Julia v1.4.2 because Atria allocates more in Julia v1.5.1."
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

    Pkg.add("PackageCompiler")
end

check_compatibility()

using PackageCompiler

cd(@__DIR__)

ver = try
    readchomp(pipeline(`grep -m1 version Project.toml`, `grep -oE '[0-9]+.[0-9]+.[0-9]+'`))
catch
    "atria"
end

if length(ARGS) == 1
    app_path = joinpath(ARGS[1], "app-$ver")
else
    app_path = joinpath(".", "app-$ver")
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
export JULIA_NUM_THREADS=8

while [[ $# -gt 0 ]]
do
    case "$1" in
        --threads)
        export JULIA_NUM_THREADS="$2"
        shift
        shift
        ;;
        -t)
        export JULIA_NUM_THREADS="$2"
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

JULIA="$DIR/julia"

"$JULIA" --color=yes -O3 --check-bounds=no --math-mode=fast "-J${DIR}/AtriaEntry.${EXT}" -e 'Atria = Base.loaded_modules[Base.PkgId(Base.UUID("226cbef3-b485-431c-85c2-d8bd8da14025"), "Atria")]; Atria.julia_main()' -- "${ALLARGS[@]}"

"""

precompile_execution_file = joinpath("test", "runtests.jl")

create_app(".", app_path, incremental=false, force=true, precompile_execution_file=precompile_execution_file, app_name="AtriaEntry")

ext = Sys.isapple() ? "dylib" : "so"
isfile(joinpath(app_path, "bin", "AtriaEntry.$ext"))

# copy julia exe
julia_exe = joinpath(Sys.BINDIR, Base.julia_exename())
cp(julia_exe, joinpath(app_path, "bin", Base.julia_exename()))

Atria_bash = joinpath(app_path, "bin", "atria")
io = open(Atria_bash, "w+")
write(io, bash_wrapper)
close(io)
chmod(Atria_bash, 0o755)

@info "Success. Atria is installed at $app_path/bin/atria"

check_compatibility()
