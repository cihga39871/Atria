
using Distributed

CJC_DEP_SAMTOOLS = get(ENV, "CJC_DEP_SAMTOOLS", "samtools")
# CJC_DEP_BEDTOOLS = get(ENV, "CJC_DEP_BEDTOOLS", "bedtools")

"""
	bam_filter(bams, output_dir::String; out_format::Symbol=:bam, flag_include_all::Int=0, flag_include_none::Int=0, flag_exclude_all::Int=0, samtools::String="$CJC_DEP_SAMTOOLS", resume::Bool=false, gzip_fastq::Bool=false, proc::Int=1)

# Parameters

- `bams::Union{String, Vector{String}}`: the path of bam file(s).

- `out_format::Symbol`: the format of output files. It should be one of `:bam`, `:sam`, `:cram`, `:fastq` (or `:fq`), `:fasta` (or `:fa`).

- `flag_include_all::Int`: only include reads with all of the FLAGs in INT present. Flag values in SAM/BAM files represent for a given combination of read pair properties. Eg: `12` means read & mate unmapped; `2` means read mapped in proper pair. (`-f` in `samtools view`)

- `flag_include_none::Int`: only include reads with none of the FLAGS in INT present. (`-F` in `samtools view`)

- `flag_exclude_all::Int`: only exclude reads with all of the FLAGs in INT present. (`-G` in `samtools view`)

- `gzip_fastq::Bool=false`: compress output `fasta` or `fastq` files. It does not work with `bam`, `sam`, `cram` files.

# Return value

- If `out_format` is one of `:bam`, `:sam`, `:cram`, return `out_files::Vector{String}`.

- If `out_format` is one of `:fastq` (or `:fq`), `:fasta` (or `:fa`): return `(out_R1s::Vector{String}, out_R2s::Vector{String})`

> **Caution**: if `out_format` is one of `:fastq` (or `:fq`), `:fasta` (or `:fa`), the input `bams` should be paired-end.

# Shell environment: modify shell environment to apply changes to the default setting:

- `CJC_DEP_SAMTOOLS`: how to call `samtools` in terminal ($CJC_DEP_SAMTOOLS)
"""
function bam_filter() end

function bam_filter(bam::String, output_dir::String; out_format::Symbol=:bam, flag_include_all::Int=0, flag_include_none::Int=0, flag_exclude_all::Int=0, samtools::String=CJC_DEP_SAMTOOLS, resume::Bool=false, gzip_fastq::Bool=false, proc::Int=1)

    ### return value are difference according to `out_format`
    if out_format in [:bam, :sam, :cram]
        empty_out = Vector{String}()
    elseif out_format in [:fastq, :fq, :fasta, :fa]
        empty_out = (Vector{String}(), Vector{String}())
    else
        throw(ArgumentError("out_format is not one of :bam, :sam, :cram, :fastq, :fq"))
        return Vector{String}()
    end

    ### validation
    try
        readchomp(`$samtools --version`)
    catch
        throw(ArgumentError("samtools not found at $samtools"))
        return empty_out
    end

    isfile(bam) || throw(ArgumentError("File $bam does not exist"))

    mkpath(output_dir)

    output_dir = realpath(output_dir)
    bam = realpath(bam)

	out_flag = (flag_include_all == 0 ? "" : "-f$flag_include_all") *
		(flag_include_none == 0 ? "" : "-F$flag_include_none") *
		(flag_exclude_all == 0 ? "" : "-G$flag_exclude_all")

    ### main code
    if out_format in [:fastq, :fq, :fasta, :fa]
		fastx = out_format in [:fastq, :fq] ? "fastq" : "fasta"
        out_r1 = basename(bam) * ".flag$out_flag.$fastx"
        out_r2 = replace(out_r1, r"_R1" => "_R2")
        if out_r1 == out_r2
            out_r1 = basename(bam) * "_R1.flag$out_flag.$fastx"
            out_r2 = basename(bam) * "_R2.flag$out_flag.$fastx"
        end

        out_r1 = joinpath(output_dir, out_r1)
        out_r2 = joinpath(output_dir, out_r2)
        complete_mark = joinpath(output_dir, ".complete-bam-filter." * basename(out_r1))

        @info(
			"Bam filter ($bam)",
			BAM_FLAG = out_flag,
			OUTPUT_R1 = out_r1 * (gzip_fastq ? ".gz" : ""),
			OUTPUT_R2 = out_r2 * (gzip_fastq ? ".gz" : "")
		)

        if resume && isfile(complete_mark) &&
            (isfile(out_r1) || isfile(out_r1*".gz")) &&
            (isfile(out_r2) || isfile(out_r2*".gz"))
            @warn "Resume: Skip completed task" TASK = "Bam filter ($bam)" BAM_FLAG = out_flag OUTPUT_R1 = out_r1 OUTPUT_R2 = out_r2
            out_r1 = isfile(out_r1) ? out_r1 : out_r1*".gz"
            out_r2 = isfile(out_r2) ? out_r2 : out_r2*".gz"

			if isfileempty(out_r1) || isfileempty(out_r2)
				@error "Bam filter ($bam): output file(s) are empty and not returned!" BAM_FLAG = out_flag OUTPUT_R1 = out_r1 OUTPUT_R2 = out_r2
				return empty_out
			else
	            return ([out_r1], [out_r2])
			end
        else
            # clear previous analysis result
            rm(out_r1, force=true)
            rm(out_r2, force=true)
            rm(out_r1*".gz", force=true)
            rm(out_r2*".gz", force=true)
            rm(complete_mark, force=true)

            try
				run(`$samtools $fastx -1 $out_r1 -2 $out_r2 -f $flag_include_all -F $flag_include_none -G $flag_exclude_all -n $bam`)
                if gzip_fastq
                    run(`gzip $out_r1 $out_r2`)
                    out_r1 *= ".gz"
                    out_r2 *= ".gz"
                end
                touch(complete_mark)

				if isfileempty(out_r1) || isfileempty(out_r2)
					@error "Bam filter ($bam): output file(s) are empty and not returned!" BAM_FLAG = out_flag OUTPUT_R1 = out_r1 OUTPUT_R2 = out_r2
					return empty_out
				else
		            return ([out_r1], [out_r2])
				end
            catch e
                rm(out_r1, force=true)
                rm(out_r2, force=true)
                rm(out_r1*".gz", force=true)
                rm(out_r2*".gz", force=true)
                rm(complete_mark, force=true)

                rethrow(e)
                @error("Bam filter failed ($bam)",
                    BAM_FLAG = out_flag,
                    OUTPUT_R1 = out_r1,
                    OUTPUT_R2 = out_r2)
                return empty_out
            end
        end
    elseif out_format in [:sam, :bam, :cram]
        out_r1 = joinpath(output_dir, basename(bam) * ".flag$out_flag.$out_format")
        complete_mark = joinpath(output_dir, ".complete-bam-filter." * basename(out_r1))

		@info(
			"Bam filter ($bam)",
			BAM_FLAG = out_flag,
			OUTPUT = out_r1
		)

        if resume && isfile(complete_mark) && isfile(out_r1)
            @warn "Resume: Skip completed task" TASK = "Bam filter ($bam)" BAM_FLAG = out_flag OUTPUT = out_r1

			if issamempty(out_r1, samtools=samtools)
				@error "Bam filter ($bam): output file is empty and not returned!" BAM_FLAG = out_flag OUTPUT = out_r1
				return empty_out
			else
	            return [out_r1]
			end
        else
            # clear previous analysis result
            rm(out_r1, force=true)
            rm(complete_mark, force=true)

            try
                run(`$samtools view --output-fmt $out_format -o $out_r1 -f $flag_include_all -F $flag_include_none -G $flag_exclude_all $bam`)
                touch(complete_mark)

				if issamempty(out_r1, samtools=samtools)
					@error "Bam filter ($bam): output file is empty and not returned!" BAM_FLAG = out_flag OUTPUT = out_r1
					return empty_out
				else
					return [out_r1]
				end
            catch e
                rm(out_r1, force=true)
                rm(complete_mark, force=true)

                rethrow(e)
                @error("Bam filter failed ($bam)",
                    BAM_FLAG = out_flag,
                    OUTPUT = out_r1)
                return empty_out
            end
        end
    end
end

function bam_filter(bams::Vector{String}, output_dir::String; out_format::Symbol=:bam, flag_include_all::Int=0, flag_include_none::Int=0, flag_exclude_all::Int=0, samtools::String=CJC_DEP_SAMTOOLS, resume::Bool=false, gzip_fastq::Bool=false, proc::Int=1)

    ### return value are difference according to `out_format`
    if out_format in [:bam, :sam, :cram]
        empty_out = Vector{String}()
    elseif out_format in [:fastq, :fq, :fasta, :fa]
        empty_out = (Vector{String}(), Vector{String}())
    else
        throw(ArgumentError("out_format is not one of :bam, :sam, :cram, :fastq, :fq"))
        return Vector{String}()
    end

    ### validation
    try
        readchomp(`$samtools --version`)
    catch
        throw(ArgumentError("samtools not found at $samtools"))
        return empty_out
    end

    foreach(f -> isfile(f) || throw(ArgumentError("File $f does not exist")), bams)

    mkpath(output_dir)

    output_dir = realpath(output_dir)

    nbam = length(bams)

    if proc == 1 || nbam == 1
        if out_format in [:bam, :sam, :cram]
            out_files = Vector{String}()
            for i in 1:nbam
                out_file = bam_filter(bams[i], output_dir; out_format=out_format, flag_include_all=flag_include_all, flag_include_none=flag_include_none, flag_exclude_all=flag_exclude_all, samtools=samtools, resume=resule, gzip_fastq=gzip_fastq)
                isempty(out_file) || append!(out_files, out_file)
            end
            return out_files
        elseif out_format in [:fastq, :fq, :fasta, :fa]
            out_r1s = Vector{String}()
            out_r2s = Vector{String}()
            for i in 1:nbam
                out_r1, out_r2 = bam_filter(bams[i], output_dir; out_format=out_format, flag_include_all=flag_include_all, flag_include_none=flag_include_none, flag_exclude_all=flag_exclude_all, samtools=samtools, resume=resume, gzip_fastq=gzip_fastq)
                isempty(out_r1) || append!(out_r1s, out_r1)
                isempty(out_r2) || append!(out_r2s, out_r2)
            end
            return (out_r1s, out_r2s)
        end
    else
        # add proc
   	   	if proc > nbam
   	   	   	proc = nbam
   	   	end

		#= I have to specify `Distributed`, otherwise compiled pcc won't work
		because of undef var error: addprocs.
		I do not know why, but it might be a bug of Distributed package, or julia itself=#
		new_proc_pids = Distributed.addprocs(proc)

        # make velvet_assembly available for every processor
        # cannot use PROGRAM_FILE!
		# @show @__FILE__
		# eval(Expr(:toplevel, :(@everywhere using Distributed)))
		if !fetch(@spawnat :any isdefined(@__MODULE__, :bam_filter))
	        @everywhere include(@__FILE__)
		end

        # pmap
   	   	output_pmap = pmap(bam -> bam_filter(bam, output_dir; out_format=out_format, flag_include_all=flag_include_all, flag_include_none=flag_include_none, flag_exclude_all=flag_exclude_all, samtools=samtools, resume=resume, gzip_fastq=gzip_fastq), bams)

   	   	rmprocs(new_proc_pids)

        # res
        if out_format in [:bam, :sam, :cram]
            out_files = Vector{String}()
            for out_file in output_pmap
                isempty(out_file) || append!(out_files, out_file)
            end
            return out_files
        elseif out_format in [:fastq, :fq, :fasta, :fa]
            out_r1s = Vector{String}()
            out_r2s = Vector{String}()
            for out_file in output_pmap
                out_r1, out_r2 = out_file
                isempty(out_r1) || append!(out_r1s, out_r1)
                isempty(out_r2) || append!(out_r2s, out_r2)
            end
            return (out_r1s, out_r2s)
        end
    end
end

function isfileempty(file::AbstractString)
	filesize(file) == 0 && return true

	# For compressed file
	# get header
    io = open(file, "r")
    head = UInt8[]
    readbytes!(io, head, 6)
    close(io)

    ###### judge compress type by file magic and then decompress
    # .gz
    if ismagicmatch(head, 1:3, UInt8[0x1f, 0x8b, 0x08])
		command = `gzip -cd`
    # .zip
    elseif ismagicmatch(head, 1:4, UInt8[0x50, 0x4b, 0x03, 0x04])
		command = `unzip -p`
    # .xz
    elseif ismagicmatch(head, 1:6, UInt8[0xfd, 0x37, 0x7a, 0x58, 0x5a, 0x00])
		command = `xz -cd`
    # .Z
    elseif ismagicmatch(head, 1:2, UInt8[0x1f, 0x9d])
		command = `zcat`
    # bzip2 .bz2
    elseif ismagicmatch(head, 1:3, UInt8[0x42, 0x5a, 0x68])
		command = `bzip2 -cd`
	else
		return false
	end

    try
		return readchomp(pipeline(`$command $file`, `head -n1`, `wc -c`)) == "0"
    catch
        @error "Decompress: command exit with error" COMMAND=`$command $file`
		# not the compress type, so the file is not empty, so return false
        return false
    end
end

function ismagicmatch(file_head::Vector{UInt8}, magic_range, expect_magic::Vector{UInt8})
    if length(file_head) >= maximum(magic_range)
        if file_head[magic_range] == expect_magic
            return true
        else
            return false
        end
    else
        @warn "Ismagicmatch: not match for `file_head` not long enough."
        return false
    end
end

function issamempty(file::String; samtools::String=CJC_DEP_SAMTOOLS)
	isfile(file) || return true

	try
		bytes_first_line = readchomp(pipeline(`$samtools view $file`, `head -n 1`, `wc -c`))
		if parse(Int, bytes_first_line) > 3
			# line not empty -> file not empty
			return false
		else
			return true
		end
	catch
		@error "issamempty() failed: command exit with error" COMMAND=`$samtools view $file`
		return false
	end
end
