
################################ 从.plr文件中提取气动数据 ######################################

# 1. 读取文本文件，提取"Polar Data"部分的内容
function read_polar_data(filepath::String)
    open(filepath, "r") do file
        data_section = false
        aoa = Float64[] # 迎角数组
        cl = Float64[]  # 升力系数数组
        cd = Float64[]  # 阻力系数数组
        cm = Float64[]  # 力矩系数数组
        
        for line in eachline(file)
            if occursin("Polar Data", line)
                data_section = true
                continue
            end
            
            # 匹配数据行并提取数据
            if data_section && occursin(r"^-?\d+\.\d+", line)
                split_line = split(line)
                if length(split_line) >= 4
                    push!(aoa, parse(Float64, split_line[1]))
                    push!(cl, parse(Float64, split_line[2]))
                    push!(cd, parse(Float64, split_line[3]))
                    push!(cm, parse(Float64, split_line[4]))
                end
            end
        end
        
        # 返回包含提取数据的DataFrame
        return DataFrame(AOA=aoa, CL=cl, CD=cd, CM=cm)
    end
end

# 2. 保存提取的数据到CSV文件
function save_to_csv(data::DataFrame, output_filepath::String)
    CSV.write(output_filepath, data)
    println("已输出文件: $output_filepath")
end

# 主函数，运行整个过程
function process_polar_data(filepath::String, output_filepath::String)
    data = read_polar_data(filepath)
    save_to_csv(data, output_filepath)
end

################################ 从.csv文件中对气动数据进行插值 ##################################

# 1. 读取数据
function load_data(filepath::String)
    return CSV.read(filepath, DataFrame)
end

# 2. 定义迎角范围
function define_aoa_range()
    aoa_range1 = -180:10:-20
    aoa_range2 = -19:1:19
    aoa_range3 = 20:10:180
    return vcat(aoa_range1, aoa_range2, aoa_range3)
end

# 3. 插值数据
function interpolate_data(data::DataFrame, aoa_range::Vector{Int})
    # 创建插值函数
    interp_cl = LinearInterpolation(data.AOA, data.CL)
    interp_cd = LinearInterpolation(data.AOA, data.CD)
    interp_cm = LinearInterpolation(data.AOA, data.CM)

    # 生成新的插值数据
    cl_vals = interp_cl.(aoa_range)
    cd_vals = interp_cd.(aoa_range)
    cm_vals = interp_cm.(aoa_range)

    # 返回新的DataFrame
    return DataFrame(AOA=aoa_range, CL=cl_vals, CD=cd_vals, CM=cm_vals)
end

# 4. 保存结果
function save_results(data::DataFrame, output_filepath::String)
    CSV.write(output_filepath, data)
    println("已输出文件: $output_filepath")
end

# 主函数运行整个过程
function process_aerodynamic_data(filepath::String, output_filepath::String)
    data = load_data(filepath)
    aoa_range = define_aoa_range()
    interpolated_data = interpolate_data(data, aoa_range)
    save_results(interpolated_data, output_filepath)
end

################################ 将插值数据转换为DUST可用的C81数据表 ##################################

# 辅助函数生成新的DataFrame并保存
# 辅助函数生成新的DataFrame并保存并提示
function process_and_save(df, col1, col2, suffix)
    # 提取相关列
    selected_df = df[:, [col1, col2]]

    # 复制第二列数据8次并插入新列
    for i in 1:8
        selected_df[:, Symbol(col2, "_copy_", i)] = selected_df[:, col2]
    end

    # 重命名第二列到第九列的列名为0.1到0.8
    for i in 1:8
        rename!(selected_df, Symbol(col2, "_copy_", i) => Symbol(string(0.1 * i)))
    end

    # 保存新的DataFrame到CSV文件
    new_file_path = replace(file_path, ".csv" => "_$suffix.csv")
    CSV.write(new_file_path, selected_df)
    println("已输出文件: $new_file_path")

    return selected_df
end

# 处理并保存CL、CD、CM数据并合并为最终文件
function process_and_return_combined(df, col1, col2, suffix)
    # 提取相关列
    selected_df = df[:, [col1, col2]]

    # 复制第二列数据8次并插入新列
    for i in 1:8
        selected_df[:, Symbol(col2, "_copy_", i)] = selected_df[:, col2]
    end

    # 重命名第二列到第九列的列名为0.1到0.8
    for i in 1:8
        rename!(selected_df, Symbol(col2, "_copy_", i) => Symbol(string(0.1 * i)))
    end

    # 统一列名并添加标识列
    rename!(selected_df, col2 => :Value)
    selected_df[!, :Type] .= suffix

    return selected_df
end

# 定义一个函数，从文件路径中提取文件名（不包含扩展名）
function extract_filename(filepath::String)
    # 使用splitpath函数分割文件路径
    path_components = splitpath(filepath)
    
    # 获取路径的最后一个元素作为文件名
    fname_with_ext = path_components[end]
    
    # 使用splitext函数去除文件扩展名
    fname_no_ext, _ = splitext(fname_with_ext)
    
    # 返回不含扩展名的文件名
    return fname_no_ext
end

# 主函数：处理输入的.plr文件，进行插值处理，并输出最终的combined_file
function dust_c81_pre(plr_filepath::String)
    # 构建文件路径
    fname = extract_filename(plr_filepath)
    output_file_path = "data/$(fname)_combined.csv"

    # 读取.plr文件并提取数据
    data = read_polar_data(plr_filepath)
    
    # 定义迎角范围
    aoa_range = define_aoa_range()

    # 对数据进行插值处理
    interpolated_data = interpolate_data(data, aoa_range)
    
    # 处理并合并CL、CD、CM数据
    cl_df = process_and_return_combined(interpolated_data, :AOA, :CL, "CL")
    cd_df = process_and_return_combined(interpolated_data, :AOA, :CD, "CD")
    cm_df = process_and_return_combined(interpolated_data, :AOA, :CM, "CM")

    # 合并CL、CD、CM数据
    combined_df = vcat(cl_df, cd_df, cm_df)
    
    # 保存合并后的数据
    CSV.write(output_file_path, combined_df)
    println("已输出文件: $output_file_path")
end