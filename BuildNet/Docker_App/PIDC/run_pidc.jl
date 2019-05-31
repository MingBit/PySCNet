using Pkg
Pkg.add("NetworkInference")
using NetworkInference

function run_pidc(Expr, filename = "pidc_links.txt")

        nodes = get_nodes(Expr)
        inferred_network = InferredNetwork(PIDCNetworkInference(), nodes)
        write_network_file(filename, inferred_network)
        end


run_pidc("100_yeast2_medium.txt")