
@inline function insert_size_decision(a_insert_size::Int64, a_score::Float64, b_insert_size::Int64, b_score::Float64; insert_size_diff::Int64 = 0)
    if a_insert_size == b_insert_size
        insert_size = a_insert_size
        score = a_score + b_score
    elseif abs(a_insert_size - b_insert_size) <= insert_size_diff
        insert_size = min(a_insert_size, b_insert_size)
        score = a_score + b_score
    elseif a_score > b_score
        insert_size = a_insert_size
        score = a_score
    else
        insert_size = b_insert_size
        score = b_score
    end
    insert_size, score
end

@inline function insert_size_decision_separate(a_insert_size::Int64, a_score::Float64, b_insert_size::Int64, b_score::Float64; insert_size_diff::Int64 = 0, score_diff::Float64 = 8.0)
    if abs(a_insert_size - b_insert_size) <= insert_size_diff
        score = a_score + b_score
        # insert sizes not changed
    elseif abs(a_score - b_score) < score_diff
        score = (a_score + b_score) / 2
        # choose the min insert size for both a and b
        if a_insert_size > b_insert_size
            a_insert_size = b_insert_size
        else
            b_insert_size = a_insert_size
        end
    elseif a_score > b_score
        score = a_score
        b_insert_size = a_insert_size
    else
        score = b_score
        a_insert_size = b_insert_size
    end
    a_insert_size, b_insert_size, score
end

@inline function is_false_positive(r1_adapter_insert_size::Int64, r1_pe_insert_size::Int64, r1_length::Int64, r2_adapter_insert_size::Int64, r2_pe_insert_size::Int64, r2_length::Int64; insert_size_diff::Int64 = 0, tail_length::Int64 = 8)::Bool

    # skip running this function when length are different.
    if r1_length < max(r2_adapter_insert_size, r2_pe_insert_size) ||
        r2_length < max(r1_adapter_insert_size, r1_pe_insert_size)
        return false
    end

    r1_adapter_error = (r1_adapter_insert_size > r1_length - tail_length) | (r1_adapter_insert_size == -1)
    r1_pe_error = r1_pe_insert_size < r1_length - tail_length
    r1_error = r1_adapter_error & r1_pe_error

    r2_adapter_error = (r2_adapter_insert_size > r2_length - tail_length) | (r2_adapter_insert_size == -1)
    r2_pe_error = r2_pe_insert_size < r2_length - tail_length
    r2_error = r2_adapter_error & r2_pe_error

    r1_adapter_inrange = abs(r1_adapter_insert_size - r1_pe_insert_size) <= insert_size_diff
    r2_adapter_inrange = abs(r2_adapter_insert_size - r2_pe_insert_size) <= insert_size_diff

    not_false_positive = r1_adapter_inrange | r2_adapter_inrange

    (r1_error | r2_error) & !not_false_positive
end
