using Pkg
Pkg.activate("InstanceGeneration")

using DataFrames, CSV

const header_str = """
[out:xml] [timeout:25];
(
"""

const end_str = """
);
(._;>;);
out body;
"""
#const bbox = " 42.43548,-71.1897,42.47413,-71.10952"
const bbox = "{{bbox}}"

function main()

    file = "InstanceGeneration/data/OSM_Mapping.csv"
    df = DataFrame(CSV.File(file))

    gdf = groupby(df,:Category)

    for (key, sub_df) in pairs(gdf)
        if key[:Category] == "N/A"
            continue
        end

        out_str = header_str
        for row in eachrow(sub_df)
            if row.Value == "*"
                out_str = out_str * "    node[\"$(row.Key)\"]($bbox);\n"
                out_str = out_str * "    way[\"$(row.Key)\"]($bbox);\n"
                out_str = out_str * "    relation[\"$(row.Key)\"]($bbox);\n"
            else
                out_str = out_str * "    node[\"$(row.Key)\"=\"$(row.Value)\"]($bbox);\n"
                out_str = out_str * "    way[\"$(row.Key)\"=\"$(row.Value)\"]($bbox);\n"
                out_str = out_str * "    relation[\"$(row.Key)\"=\"$(row.Value)\"]($bbox);\n"
            end
        end

        out_str = out_str * end_str
        println("Start Query Script for $key\n")
        println(out_str)
        println("\nEnd Query Script for $key")
    end

end

main()