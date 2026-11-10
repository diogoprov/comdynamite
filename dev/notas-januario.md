# Poder do teste de Januario et al. (2025)

Januario et al. (2025, *Ecology Letters*, doi:10.1111/ele.70021, conferido no texto completo,
PMC11687345) não acharam relação entre a persistência metapopulacional (ocupação de equilíbrio ε,
resiliência τ) de 189 peixes demersais e a taxa de especiação nas pontas (λ_DR, BAMM), por PGLS.
Pergunta: se a persistência causasse a especiação, esse desenho detectaria?

Scripts: `data-raw/januario_power.R` (rodada principal, `januario_power.rds`) e
`data-raw/januario_power_ols.R` (diagnóstico). Desenho em cada réplica: trait latente X por BM,
árvore de nascimento–morte com λ = exp(bX), 200 pontas; X padronizado define a dispersão de cada
espécie em `simulate_metacommunity()` (grade 10 × 10, sem competição); ε e τ estimados por um
modelo markoviano de duas transições nos últimos `window` passos; PGLS (`nlme::gls` +
`ape::corBrownian`) de log λ_DR, e do log λ verdadeiro, contra ε + τ.

## Rodada principal (200 réplicas por célula; 44 min em 6 núcleos; 0 erros)

| b | janela | cor(τ, z) | cor(ε, z) | cor(λ_DR, λ) | rejeição, PGLS λ_DR | rejeição, PGLS λ verdadeiro |
|---|---|---|---|---|---|---|
| 0 | 10 / 30 / 100 | 0,86–0,88 | 0,34–0,37 | — | 0,09 / 0,08 / 0,04 | — |
| 0,25 | 10 / 30 / 100 | 0,85–0,87 | 0,34–0,37 | 0,33–0,35 | 0,085 / 0,09 / 0,08 | 0,97 / 1,00 / 1,00 |
| 0,5 | 10 / 30 / 100 | 0,84–0,87 | 0,33–0,36 | 0,50–0,52 | 0,08 / 0,075 / 0,08 | 0,97 / 0,99 / 1,00 |
| 1 | 10 / 30 / 100 | 0,83–0,86 | 0,33–0,35 | 0,61–0,64 | 0,075 / 0,06 / 0,05 | 1,00 / 0,99 / 0,99 |

(correlações de Spearman, medianas entre réplicas; teste conjunto ε + τ a 5%)

Controle negativo (X → extirpação, janela 30, 100 réplicas): ε e τ não acompanham z
(cor ≈ −0,06 e 0,03–0,08), e a PGLS com λ verdadeiro rejeita 1% com b = 1; com λ_DR, 13% com
b = 0 e 4% com b = 1.

## Diagnóstico (60 réplicas por b, dispersão, janela 30)

| b | PGLS | OLS | Spearman λ_DR × τ | Spearman λ_DR × z |
|---|---|---|---|---|
| 0 | 0,05 | 0,20 | 0,15 | 0,20 |
| 1 | 0,05 | 1,00 | 1,00 | 1,00 |

## Leitura

- ε e τ são bem estimados (τ melhor que ε) mesmo com janela de 10 passos: o ruído da persistência
  não é o gargalo.
- A PGLS com λ_DR não tem poder nenhum, em nenhuma combinação: a taxa de rejeição fica no nível do
  tamanho (5–9%) mesmo com b = 1, quando λ_DR tem correlação 0,6 com o λ verdadeiro. A mesma PGLS
  com o λ verdadeiro detecta em ~100%.
- Sem correção filogenética, λ_DR detecta a ligação em 100% das réplicas, mas com tamanho inflado
  (15–20% sob b = 0).
- Interpretação nossa (não testada à parte): λ_DR é uma média ponderada ao longo do caminho até a
  raiz, então espécies irmãs têm λ_DR quase iguais e o sinal fica todo entre clados. A PGLS
  browniana trabalha com contrastes, que são justamente onde λ_DR não carrega informação. Com X
  evoluindo por BM, a associação X–λ_DR também é entre clados. O resultado: o teste que controla
  o tamanho não vê nada, e o teste que vê não controla o tamanho.
- Implicação para Januario et al.: a ausência de relação na PGLS com λ_DR é o resultado esperado
  sob esse modelo mesmo com ligação causal forte. O "desacoplamento" não distingue ausência de
  ligação de falta de poder. Testes de correlação com taxas de ponta com nulo por simulação
  (ES-sim; Harvey & Rabosky 2018, *Methods Ecol. Evol.* 9:984–993,
  doi:10.1111/2041-210X.12949, conferido no Crossref) são a alternativa natural; testado abaixo.
- Ressalvas: a ligação X → λ é imposta (SSE contínuo); dispersão como via de persistência; 200
  espécies sem filtro de informação; sem detecção imperfeita; BAMM não testado.

## ES-sim no mesmo desenho (`data-raw/januario_essim.R`, `januario_essim.rds`)

ES-sim implementado da descrição de Harvey & Rabosky (2018): correlação de Pearson entre log
λ_DR e o trait, nulo por traits brownianos na árvore (L z; a correlação não muda com a raiz e a
taxa do BM), p bilateral com 999 sorteios. Janela 30, 200 réplicas por b (dispersão) e 100 no
controle negativo (extirpação); 8,8 min em 6 núcleos, 0 erros. Mesmas réplicas para PGLS e
ES-sim.

| ligação | b | PGLS λ_DR (ε + τ) | ES-sim τ | ES-sim ε | ES-sim τ ou ε (Bonferroni) |
|---|---|---|---|---|---|
| dispersão | 0 | 0,095 | 0,015 | 0,005 | 0,01 |
| dispersão | 0,25 | 0,065 | 0,585 | 0,015 | 0,41 |
| dispersão | 0,5 | 0,065 | 0,89 | 0,045 | 0,77 |
| dispersão | 1 | 0,055 | 0,975 | 0,035 | 0,92 |
| extirpação | 0 | 0,09 | 0,01 | 0 | 0 |
| extirpação | 1 | 0,02 | 0 | 0 | 0 |

- O ES-sim com τ detecta a ligação que a PGLS não vê (0,59–0,98) e controla o tamanho (0,015,
  conservador); o controle negativo não rejeita.
- Com ε não há poder: ε acompanha mal o trait latente (Spearman 0,33–0,37).
- Implicação: no desenho de Januario et al. (2025), um teste por simulação sobre τ teria
  detectado uma ligação causal moderada; o resultado nulo da PGLS não é evidência de
  desacoplamento. Ressalvas acima continuam (ligação imposta, dispersão como via, sem detecção
  imperfeita, BAMM não testado).
