function best_linear_code_in_codetables(n::Int, k::Int, F::String; verbose::Bool=false)
    q_match = match(r"GF\((\d+)\)", F)
    if q_match === nothing
        throw(ArgumentError("Invalid finite field format. Expected format: GF(q), where q is an integer."))
    end
    q = parse(Int, q_match.captures[1])
    if !(q in [2, 3, 4, 5, 7, 8, 9])
        throw(ArgumentError("q (= $q) must be in [2, 3, 4, 5, 7, 8, 9]"))
    end
    base_url = "http://www.codetables.de/BKLC/BKLC.php"
    params = "?q=$q&n=$n&k=$k"
    url = base_url * params
    if verbose
        println("Querying URL: $url")
    end
    response = HTTP.get(url)
    if response.status != 200
        throw(IOError("Failed to retrieve data: HTTP status $(response.status)"))
    end
    html_content = String(response.body)
    document = parsehtml(html_content)
    pre_tags = eachmatch(sel"pre", document.root)
    if isempty(pre_tags)
        throw(IOError("Error parsing data: <PRE> tags not found in the response."))
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
