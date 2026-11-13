# Simulador de órbitas

![CI](https://github.com/RodrigoJBotelho/simulador-orbitas/actions/workflows/ci.yml/badge.svg)

Problema de dois corpos (satélite em torno da Terra) integrado com quatro métodos numéricos e comparado com a solução exata de Kepler. Mostra, na prática, por que métodos simpléticos conservam melhor a energia em simulações longas.

| Método | Ordem | Comportamento na energia |
|---|---|---|
| `euler` (explícito) | 1 | cresce sem limite: a órbita "abre em espiral" |
| `simpletico` (Euler simplético) | 1 | oscila, mas fica limitada |
| `verlet` (Velocity Verlet) | 2 | oscila com amplitude pequena, limitada |
| `rk4` (Runge–Kutta clássico) | 4 | erro muito pequeno, mas com deriva lenta em tempos longos |

> Projeto didático, em plano orbital 2D e com a Terra como ponto de massa. Não inclui achatamento (J2), arrasto, Lua nem Sol.

## Como compilar e rodar

Precisa apenas de um compilador C++11 (`g++` ou `clang++`) e do `make`.

```bash
make run        # compila e gera os CSVs na pasta saida/
make teste      # verifica ordens de convergência e conservação de energia
make exemplos   # copia os CSVs para docs/exemplos/ (usados pela página web)
```

Sem `make`:

```bash
g++ -O2 -std=c++11 -o orbitas orbitas/orbitas.cpp
./orbitas 0.5 5 200      # excentricidade, nº de órbitas, passos por órbita
./orbitas --verificar    # testes de sanidade (0,1 <= e <= 0,8)
```

Órbita padrão: semieixo maior de 24 000 km, excentricidade 0,5 (perigeu 12 000 km, apogeu 36 000 km, período de cerca de 10,3 h), partindo do periapsis.

## Saídas

| Arquivo | Conteúdo |
|---|---|
| `orbita_trajetoria.csv` | posições (km) de cada método e da solução exata (`x_*`, `y_*`) |
| `orbita_energia.csv` | erro relativo de energia × tempo (em órbitas) |
| `orbita_convergencia.csv` | erro de posição após 1 órbita × passos por órbita |
| `orbita_convergencia_energia.csv` | erro máximo de energia em 1 órbita × passos por órbita |

## Visualização

Abra `docs/index.html` no navegador e use **Abrir meu CSV**. Colunas `x_nome`/`y_nome` viram trajetórias (com a Terra desenhada no foco); os outros arquivos viram curvas, com escala logarítmica opcional. Ativando o GitHub Pages (Settings → Pages → pasta `/docs`), os botões de exemplo funcionam online.

## Resultados com os parâmetros padrão (5 órbitas, 200 passos/órbita)

| Método | Erro máx. de energia | Erro de posição final |
|---|---|---|
| euler | 7,4 × 10⁻¹ | ~173 000 km |
| simpletico | 5,2 × 10⁻² | ~12 300 km |
| verlet | 2,7 × 10⁻³ | ~5 300 km |
| rk4 | 5,3 × 10⁻⁶ | ~6,7 km |

Ordens observadas (erro de energia, 800 → 1600 passos/órbita): euler ≈ 0,8, simpletico ≈ 1,0, verlet ≈ 2,0, rk4 ≈ 4,1.

Observação: a ordem aparente depende da grandeza medida. No erro de **posição** após uma órbita completa (partindo do periapsis), o Euler simplético mostra ordem próxima de 2, mas a ordem do método (vista no erro de energia) é 1. Por isso o teste automático usa a energia.

## Testes automáticos

`make teste` (e o workflow em `.github/workflows/ci.yml`, que roda a cada `push`) confere:
1. as ordens de convergência dos quatro métodos;
2. que a solução exata de Kepler conserva a energia (< 10⁻⁹);
3. que o Euler explícito é pior que o Verlet na energia.

O teste é definido para 0,1 ≤ e ≤ 0,8: em órbita circular e em órbitas muito excêntricas as ordens aparentes fogem do padrão.

## Próximos passos possíveis

Passo adaptativo (Dormand–Prince), perturbação J2, três corpos restrito (Terra–Lua–satélite), integração direta com a mesma saída VTK do solver de Poisson para ver no ParaView.

## Estrutura

```
orbitas/orbitas.cpp
docs/index.html        # visualizador (HTML + JS puro)
docs/exemplos/         # CSVs de exemplo
docs/teoria.md
.github/workflows/ci.yml
Makefile  LICENSE  README.md
```

## Referências

- Hairer, E.; Lubich, C.; Wanner, G. *Geometric Numerical Integration*. Springer.
- Curtis, H. D. *Orbital Mechanics for Engineering Students*. Butterworth-Heinemann.
- Press, W. H. et al. *Numerical Recipes*. Cambridge University Press.

## Licença

MIT. Veja `LICENSE`.
