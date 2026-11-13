# Resumo da teoria

## Equação do movimento
`d²r/dt² = −GM r / |r|³`. A energia específica `ε = v²/2 − GM/r` e o momento angular `L = x·vy − y·vx` são constantes do movimento.

## Solução exata (Kepler)
Com semieixo `a`, excentricidade `e` e movimento médio `n = √(GM/a³)`:
1. anomalia média `M = n·t`;
2. equação de Kepler `E − e·sin E = M` (resolvida por Newton);
3. `x = a(cos E − e)`, `y = a√(1−e²)·sin E`.

## Métodos
- **Euler explícito:** `x⁺ = x + Δt·v`, `v⁺ = v + Δt·a(x)`. A energia aumenta a cada passo.
- **Euler simplético:** atualiza `v` primeiro e usa o `v` novo em `x`. Preserva a estrutura simplética; a energia oscila de forma limitada.
- **Velocity Verlet:** `x⁺ = x + Δt·v + ½Δt²·a`, `v⁺ = v + ½Δt(a + a⁺)`. Simplético, reversível no tempo, ordem 2.
- **RK4:** quatro avaliações da aceleração por passo, ordem 4. Não é simplético: em tempos muito longos a energia deriva lentamente.

## Convergência
Se o erro global é `C·Δtᵖ`, a ordem observada entre dois refinamentos é `p ≈ log₂(erro(N) / erro(2N))`.
