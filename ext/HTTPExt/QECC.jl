function query_quantum_code_in_codetables(n::Int, k::Int, q::Int; verbose::Bool=false)
    if !(q in [2, 3, 4, 5, 7, 8])
        throw(ArgumentError("q (= $q) must be in [2, 3, 4, 5, 7, 8]"))
    end
    if !(1 ≤ n ≤ 256)
        throw(ArgumentError("n (= $n) must be in the range 1 ≤ n ≤ 256"))
    end
    if !(0 ≤ k ≤ n)
        throw(ArgumentError("k (= $k) must be in the range 0 ≤ k ≤ n"))
    end

    # Correct base URL for quantum codes
    base_url = "https://codetables.de/QECC.php"
    params = "?q=$q&n=$n&k=$k"
    url = base_url * params

    if verbose
        println("Querying URL: $url")
    end

    response = HTTP.get(url)
    if response.status != 200
        error("Failed to retrieve data: HTTP status $(response.status)")
    end

    html_content = String(response.body)
    document = parsehtml(html_content)
    pre_tags = eachmatch(sel"pre", document.root)
    if isempty(pre_tags)
        error("Error parsing data: <PRE> tags not found in the response.")
    end

    pre_content = nodeText(pre_tags[1])
    return pre_content
end

function nodeText(node)
    text = ""
    for child in children(node)
        if isa(child, HTMLText)
            text *= child.text
        else
            text *= nodeText(child)
        end
    end
    return text
end
