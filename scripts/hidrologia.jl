using CSV, DataFrames, Plots, Dates, Statistics, BenchmarkTools

dados = DataFrame()

files = readdir("dados/")

for file ∈ files
    append!(dados, CSV.read(joinpath("dados", file), DataFrame; select = [
        :nom_reservatorio,
        :din_instante,
        :val_vazaoafluente,
        :val_vazaoturbinada
    ]))
end

usinas = [
    "SOBRAGI",
    "TRÊS MARIAS"
]

function transform_bymonth_mean(dados::DataFrame, usina::String)
    dados |>
    df -> subset(df, :nom_reservatorio => x -> x .== usina) |>
    df -> transform(df, :din_instante => ByRow(yearmonth) => :din_instante) |>
    df -> groupby(df, :din_instante)
end

function afluencia_mean(dados::DataFrame, usina::String)
    dados |>
    df -> transform_bymonth_mean(df, usina) |>
    df -> combine(df, :val_vazaoafluente => mean, :val_vazaoturbinada => mean)
end 

sobragi_afluencia = afluencia_mean(dados, "SOBRAGI")

function get_afluencia(dados::DataFrame, usinas::Array{String}, col=:val_vazaoafluente)
    ret = DataFrame()

    for usina ∈ usinas
        s = Symbol(usina)
        
        column = dados |>
        df -> subset(df, :nom_reservatorio => x -> x .== usina) |>
        df -> transform!(df, :din_instante => ByRow(yearmonth) => :din_instante) |>
        df -> groupby(df, :din_instante) |>
        df -> combine(df, col => mean => s) |>
        df -> select!(df, s)
        
        ret = hcat(ret, column)
    end

    ret
end

afluencia = get_afluencia(dados, ["SOBRAGI", "SÃO SIMÃO", "TRÊS MARIAS"])