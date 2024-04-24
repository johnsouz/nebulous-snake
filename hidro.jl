
using DelimitedFiles, Printf
using JuMP, HiGHS, Plots

function maplinear(x, x1, y1, x2, y2)
    (x - x1) * (y2 - y1) / (x2 - x1) + y1
end

function get_hidro_data()
    afluencia_matrix = readdlm("db/afluencia.csv", ',')
    hidro_matrix = readdlm("db/hidro.csv", ',')  

    rows, cols = size(afluencia_matrix)
    n_meses = rows-1
    
    hidro_set = Symbol.(afluencia_matrix[1,:])
    colunas = Symbol.(hidro_matrix[1, 2:end])
    
   (    hidro_set = hidro_set,
        n_meses = n_meses,
        afluencia = Containers.DenseAxisArray(
            Matrix{Float64}(afluencia_matrix[2:end,:]),
            1:n_meses,
            hidro_set),
        hidro = Containers.DenseAxisArray(
            Matrix{Float64}(hidro_matrix[2:end, 2:end]),
            hidro_set,
            colunas))
end

function get_termo_data()
    termo_matrix = readdlm("db/termos.csv", ',')

    Symbol.(termo_matrix[2:end, 1]),
    Containers.DenseAxisArray(
        Matrix{Float64}(termo_matrix[2:end,2:end]),
        Symbol.(termo_matrix[2:end,1]),
        Symbol.(termo_matrix[1,2:end]))
end

hidro_set, n, afluencia, hidro_data = get_hidro_data()
termo_set, termo_data = get_termo_data()

demanda = 2200.0
descontinua = false
μ_chuva = 0.95

modelo = Model(HiGHS.Optimizer)
# modelo = Model(Clarabel.Optimizer)

# region Hidroeletricas

if descontinua
    @variable(modelo, turbinado[i ∈ 1:n, ug ∈ hidro_set] ∈ Semicontinuous(hidro_data[ug, :turbinado_min], hidro_data[ug, :turbinado_max]))
else
    @variable(modelo, hidro_data[ug, :turbinado_min] <= turbinado[i ∈ 1:n, ug ∈ hidro_set] <= hidro_data[ug, :turbinado_max])
end

@variable(modelo, 0 <= vertido[i ∈ 1:n, ug ∈ hidro_set])
@variable(modelo, hidro_data[ug, :reservatorio_min] <= armazenado[i ∈ 1:n, ug ∈ hidro_set] <= hidro_data[ug, :reservatorio_max])
@variable(modelo, 0 <= deficit[i ∈ 1:n])

#endregion
# region Termoeletricas

@variable(modelo, termo_data[ug, :geracao_min] <= geracao_termica[i ∈ 1:n, ug ∈ termo_set] <= termo_data[ug, :geracao_max])

#endregion


@constraint(modelo, balanco_hidrico_inicial[ug ∈ hidro_set],
    +    vertido[1, ug]
    + armazenado[1, ug]
    +  turbinado[1, ug]
    ==
    + 0.5 * (hidro_data[ug, :reservatorio_max] - hidro_data[ug, :reservatorio_min]) + hidro_data[ug, :reservatorio_min]
    + afluencia[1, ug] * μ_chuva
)

@constraint(modelo, balanco_hidrico[i ∈ 2:n, ug ∈ hidro_set],
    +    vertido[i, ug]
    + armazenado[i, ug]
    +  turbinado[i, ug]
    ==
    + armazenado[i-1, ug]
    +  afluencia[i, ug] * μ_chuva
)

@constraint(modelo, balanco_hidrico_final[ug ∈ hidro_set],
    armazenado[n, ug] >= .4 * (hidro_data[ug, :reservatorio_max] - hidro_data[ug, :reservatorio_min]) + hidro_data[ug, :reservatorio_min])

@constraint(modelo, balanco_energetico[i ∈ 1:n],
    + sum(turbinado[i, ug]*μ for ug ∈ hidro_set, μ ∈ hidro_data[ug, :conversao_hidrica])
    + sum(geracao_termica[i, ug] for ug ∈ termo_set)
    + deficit[i] == demanda
)

@expression(modelo, penalizacao_reservatorio[i ∈ 1:n, ug ∈ hidro_set],
    maplinear(armazenado[i, ug], hidro_data[ug, :reservatorio_max], 0, hidro_data[ug, :reservatorio_min], 1)
)

@objective(modelo, Min,
    + 1e6 * sum(deficit)
    + 1e6 * sum(vertido)
    + sum(penalizacao_reservatorio)
    + sum(geracao_termica[i, ug]*termo_data[ug, :custo] for i ∈ 1:n, ug ∈ termo_set))

optimize!(modelo)
solution_summary(modelo; verbose=true)

function plot_armazenamento(armazenado, hidro)
    p = plot(title="Volume útil do reservatório", xlabel="Mês", ylabel="%", ylims=[0.0, 1.0])
    for ug ∈ axes(hidro)[1]
        series = maplinear.(value.(armazenado[:, ug]), hidro[ug, :reservatorio_min], 0, hidro[ug, :reservatorio_max], 1)
        plot!(p, series.data, labels=String(ug))
    end
    p
end

function plot_turbinado(turbinado, hidro, percent=false)
    p = plot(xlabel="Mês", ylabel="m³⋅mês⁻¹", title="Turbinamento")
    for ug ∈ axes(hidro)[1]
        series = percent ?
            maplinear.(value.(turbinado[:, ug]), 0, 0, hidro[ug, :turbinado_max], 1) :
            value.(turbinado[:, ug])

        plot!(p, series.data;
            label = String(ug),
            legend = :outertop,
            legend_column = -1)
    end
    p
end

function plot_afluencia(afluencia)
    p = plot(xlabel="Mês", ylabel="m³⋅mês⁻¹", title="Afluencia")
    for ug ∈ axes(afluencia)[2]
        plot!(p, afluencia[:, ug].data;
            label =String(ug),
            legend = :topleft)
    end
    p
end

function plot_defluencia_filled(afluencia, turbinado, vertido, ug)
    p = plot(title=@sprintf("Balanço hídrico para ug %s", ug),
        xlabel="Mês",
        ylabel="m³⋅mês⁻¹")
    
    # afluencia
    areaplot!(p, afluencia[:, ug].data;
        label = "Afluencia",
        fillcolor = :lightskyblue,
        falpha = 0.25)
    
    # defluencia
    areaplot!(p, -hcat(value.(vertido[:, ug]).data, value.(turbinado[:, ug]).data);
        labels = ["Vertido" "Turbinado"],
        fillrange = 10,
        seriescolor = [:red :darkseagreen1],
        falpha = 0.25)
    
    plot!(p, afluencia[:, ug].data - value.(turbinado[:, ug]).data - value.(vertido[:, ug]).data;
        label = "Diferença",
        color = :black,
        linewidth = 1.5)
    p
end


function plot_defluencia_lines(afluencia, turbinado, vertido, ug = :A)
    p = plot(xlabel="Mês", ylabel="m³⋅mês⁻¹", title="Balanço hídrico")
    
    afluencia_val = afluencia[:, ug].data
    turbinado_val = value.(turbinado[:, ug]).data
    vertido_val = value.(vertido[:, ug]).data

    # afluencia
    plot!(p, afluencia_val;
    label = "Afluente",
    color = :lightskyblue)
    
    plot!(p, -(vertido_val + turbinado_val);
    label = "Vertido",
    color = :lightsalmon)
    
    plot!(p, -turbinado_val;
    label = "Turbinado",
    color = :darkseagreen1)
    
    plot!(p, afluencia_val - turbinado_val - vertido_val)
    p
end

function plot_geracao_tipos(turbinado, termica, deficit)
    p = plot(title="Geração por tipo de fonte")
    serie_hidro = value.(turbinado)
    serie_termo = value.(termica)

    plot!(p, collect(sum(serie_hidro[i,:].data) for i ∈ axes(serie_hidro)[1]),
        label="Hidro")
    
    plot!(p, collect(sum(serie_termo[i,:].data) for i ∈ axes(serie_termo)[1]),
        label="Termo")
    
    plot!(p, value.(deficit),
        label="Deficit")

end

function plot_model()
    savefig(plot_afluencia(afluencia), "plots/afluencia.png")

    for u in hidro_set
        p = plot_defluencia_filled(afluencia, turbinado, vertido, u)
        savefig(p, "plots/defluencia_" * string(u) * ".png")
    end

    savefig(plot_turbinado(turbinado, hidro_data), "plots/turbinado.png")
    savefig(plot_armazenamento(armazenado, hidro_data), "plots/armazenamento.png")
    savefig(plot_geracao_tipos(turbinado, geracao_termica, deficit), "plots/geracao.png")   
end


# savefig(plot_afluencia(afluencia), "plots/Afluencia.png")
