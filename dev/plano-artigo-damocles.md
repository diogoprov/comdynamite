# Plano: artigo sobre montagem local e regimes de diversificação (DAMOCLES)

Escrito em 2026-11-09; atualizado no mesmo dia com as decisões do Diogo e a etapa E1
(seção 4). Documento interno de trabalho. Cada bloco de simulação e de análise empírica só
começa depois de aprovado. Não confundir com
`dev/notas-artigo-taxas.md`, que trata de taxas de evolução de traits.

## 1. Pergunta e enquadramento

**Pergunta.** As taxas de colonização (γ) e de extinção local (μ) das comunidades locais
diferem entre regimes de diversificação, além do que a própria história de especiação já
explica?

**Por que agora.** Reijenga et al. (2026) estenderam o DAMOCLES para vários clados ao mesmo
tempo. Eles mostraram que clados de diversificação rápida têm proporcionalmente menos
espécies em simpatria e atribuíram isso só à história de especiação: γ e μ eram constantes
entre clados, e o modelo global venceu o modelo por família por AIC. Ninguém testou se γ e μ
mudam *dentro* de uma árvore, acompanhando regimes ou taxas de diversificação estimados por
outro método. É essa a ponte que o pacote passou a fazer (PRs #39, #40, #44 e #56), com:

- regimes e covariáveis de ramo em `damocles_fit()`;
- ajustes sobre configurações do posterior em `damocles_samples()`;
- teste contra partições aleatórias em `damocles_partition_test()`;
- taxas de ponta e de ramo com `dr_rates()` e os exemplos de BAMM e ClaDS.

A busca bibliográfica de 2026-11-09 (OpenAlex, Scholar Gateway e Undermind) não encontrou
precedente para essa integração.

**Precedente por passo:**

- o ajuste compartilhado já tem precedente direto (Reijenga et al. 2026);
- a simpatria tem análogos em DEC e GeoSSE;
- LRT e intervalos por perfil são inferência padrão;
- a novidade está em regimes e covariáveis vindos de BAMM e ClaDS, na propagação do
  posterior, no teste de partições e na diluição por regressão.

**Formato (decidido em 2026-11-09): artigo único para a AmNat.** Os métodos (seção 3)
entram como suporte, com as simulações em apêndice, e a pergunta biológica conduz o texto.
Consequência: Benício et al. (2021) sozinho, com n = 14, não sustenta a inferência de um
artigo da AmNat. A escala da Mata Atlântica (E6, com o ATLANTIC AMPHIBIANS) passa a fazer
parte do artigo, e Benício vira o caso com amostragem padronizada e contraste costa contra
interior. O risco de revisores de métodos e de ecologia cobrarem coisas diferentes continua;
a mitigação é um apêndice de simulações completo e um texto principal centrado na pergunta.

## 2. Dados

Tudo em `data-raw/` no checkout principal, sem commit. A pasta `data-raw` está no
`.Rbuildignore`. Licenças a respeitar antes de qualquer redistribuição.

| Conjunto | Onde | Conteúdo | Situação |
|---|---|---|---|
| Benício et al. (2021) | `data-raw/benicio2021/` (Dryad, doi:10.5061/dryad.sf7m0cg65) | presença e ausência de anuros em 14 áreas protegidas (7 no interior, 7 na costa), coordenadas, `geo_frogs.txt`, modos reprodutivos, scripts (GeoSSE, BioGeoBEARS) | baixado |
| Árvores de Benício | mesma pasta | `tree_pyron_112_species` e `tree_feng_112_species` (casam 112/112 com a matriz); `tree_pyron_248_species` e `tree_feng_248_species` (NEXUS, 97 árvores), com 95–96/112 por sinônimos (*Hypsiboas* → *Boana*, *Scinax*) | baixado |
| Portik et al. (2023) | `data-raw/portik2023/` (OSF 3cbsh) | árvore datada por TreePL, 5.326 pontas, binária e ultramétrica, mais 100 bootstraps datados | baixado; o OSF só diz "copyright Daniel Portik, 2023", sem licença aberta: **não redistribuir** |
| Vasquez-Restrepo (2024, preprint) | Figshare doi:10.6084/m9.figshare.24894363 (CC BY 4.0) | BAMM em 100 árvores totalmente amostradas de anfíbios; DR; presença e ausência da IUCN | só a parte de anfíbios, extraída do zip de 6,7 GB por intervalos de bytes, em `data-raw/vasquez2024/` (cerca de 1,2 GB): event data, `mcmc_out` e árvore de cada uma das 100 corridas, `Amphibians_SAMPLED TREES.nex`, taxas por espécie (BAMM, ES, FP) e comprimentos de ramo |
| ATLANTIC AMPHIBIANS (Vancine et al. 2018) | `af-frogs-climate-rarity/data/raw/atlantic_amphibians/`; também no `atlanticr` | 1.163 sítios, 17.626 registros | disponível; esforço amostral heterogêneo |

### 2.1 Casamento de nomes com a árvore do Portik (feito em 2026-11-09)

- Feito por `data-raw/benicio2021_trees.R`, que escreve `data-raw/benicio2021_taxonomy.csv`
  (nome na árvore de 248, nome na matriz, nome canônico, nome no Portik, ação, nota). O nome
  canônico é o do Portik quando existe, senão um nome corrigido; a matriz casa com as árvores
  por `comm_name` → `name`.
- Das 248: 202 casam direto, 34 por renomeação (gênero gramatical em *Boana*, grafias,
  *Hylomantis granulosa*, *Phantasmarana boticariana*, *Scinax perereca*) e 12 são
  enxertadas, 10 delas nas comunidades: *Adenomera* sp. J, *Aplastodiscus* sp. 5,
  *Bokermannohyla izecksohni*, *Dendrophryniscus* sp., *Fritziana* aff. *fissilis*,
  *Ischnocnema* sp1–sp3, *Phantasmarana* aff. *boticariana*, *Paratelmatobius* aff.
  *cardosoi*, *Phrynomedusa dryade*, *Rhinella fernandezae*.
- Enxerto: no ramo-tronco do ancestral comum (no Portik) das espécies do pool que formam o
  clado-irmão na árvore de Benício, na mesma altura relativa entre esse ancestral e o pai.
  Sempre sobre o pool de 248, que contém os clados-irmãos de todas, e depois podado.
- *Ischnocnema guentheri* C1, C3 e C4 (decidido em 2026-11-09): linhagens ainda sem nome e
  com divergências antigas, tratadas como *Ischnocnema* sp1, sp2 e sp3. Enxertadas como
  táxons distintos; separam-se de *I. henselii* a cerca de 20–23 Ma na árvore resultante.

### 2.2 Definição do pool regional (decisão central)

- Com a árvore de 112 espécies, o pool é a união das 14 comunidades: nenhuma espécie está
  ausente de todos os sítios. Isso condiciona o DAMOCLES a linhagens que colonizaram algum
  lugar e tende a superestimar γ.
- **O pool de 248 não é um pool regional.** O README do Dryad descreve o alinhamento como
  as sequências das 248 espécies usadas no RAxML, e o Diogo suspeita que as 136 espécies a
  mais entraram por terem sequências no GenBank, para apoiar a inferência. A lista confirma:
  inclui espécies do Nordeste, como *Adelophryne baturitensis*, e da Amazônia. Serve como
  análise de sensibilidade, não como pool principal.
- **Pool regional (decidido em 2026-11-09): Mata Atlântica da Serra do Mar e do interior,
  incluindo as regiões de Cerrado.** Operacionalizado como os anuros com nome válido nos 447
  sítios do ATLANTIC AMPHIBIANS no estado de São Paulo, que cobrem costa, Serra do Mar e
  interior, e contêm as 14 áreas. Ficam de fora *Rana catesbeiana* (introduzida) e
  *I. guentheri*: os registros paulistas do complexo são tomados como as linhagens sp1–sp3,
  que entram pelas comunidades, para não contar o complexo duas vezes. Somam-se as 10
  espécies das comunidades sem nome válido no ATLANTIC AMPHIBIANS. Total: 222 espécies
  (`data-raw/benicio2021_regional_pool.csv`).
  - 199 casam com o Portik (com *Ololygon* → *Scinax*, *Pithecopus* → *Phyllomedusa*,
    *Phantasmarana* e gênero do epíteto);
  - 9 são enxertos guiados pela árvore de Benício;
  - 14 são enxertadas em ponto uniforme sobre os ramos do clado dos congêneres, sorteado de
    novo em cada árvore. Para as ex-*Ololygon*, o clado é o grupo *catharinae* (21 pontas,
    32 Ma). A mais frequente é *Scinax hiemalis* (41 sítios); as demais têm de 1 a 6 sítios.
  - Limite: os registros sem nome válido (526 em SP) ficam de fora; "sp." e "aff." de estudos
    diferentes não são comparáveis.
  - Falsas ausências: o Diogo lembrou que o ATLANTIC AMPHIBIANS seleciona estudos para
    reduzir falsas ausências. O que os metadados (Vancine et al. 2018, MetadataS1, seção
    "Research methods") descrevem como critério de inclusão:
    - entram só estudos que amostraram comunidades de anfíbios;
    - ficam fora descrições de espécies e trabalhos de distribuição de um único táxon;
    - ficam fora estudos sem composição de espécies detalhada ou sem coordenadas.

    Não encontrei nos metadados limiar de duração nem de método. O esforço é registrado
    como a duração total do estudo em meses: média de 14,3 meses, desvio de 14,1, de menos
    de 1 a 191. Os autores dizem que essa medida superestima o esforço real, e a relação
    entre riqueza e esforço é fraca (R² = 0,12). Os campos `effort_months`, `active_methods`
    e `passive_methods` permitem uma análise de sensibilidade com filtros mais estritos.
- Plano:
  - pool regional como análise principal;
  - pools de 112 e de 248 como sensibilidade;
  - a simulação S6 quantifica o viés esperado.
- Costa e interior: o DAMOCLES assume um único pool. Duas formas, a comparar:
  - (a) pool combinado com região como regime de γ e μ nos ramos terminais, o que exige
    estender `regimes` para depender da comunidade, não só do ramo (ver 3.8);
  - (b) ajustes separados por região, comparados por LRT contra o modelo compartilhado.

## 3. Lacunas metodológicas e simulações

Todas em `data-raw/`, um script por bloco, com:

- `RNGkind("L'Ecuyer-CMRG")`, `CORES` por variável de ambiente e `.rds` por cenário
  (pontos de retomada);
- `message()` de progresso por cenário;
- nenhum `repeat` sem limite de tentativas (lição da seção 7 da validação, que travou por
  cerca de 6 h).

Árvores simuladas com `treeshapesim` ou com processos de nascimento e morte com mudanças de
taxa (pacotes GPL só em `data-raw/`).

Os tempos abaixo são **estimativas grosseiras**, a medir num piloto antes de cada bloco.

### S1. Erro tipo I com diversificação heterogênea (prioridade máxima)

- **Pergunta:** o teste de regime rejeita γ constante só porque a árvore tem mudanças de
  taxa de especiação?
- **Desenho:**
  - árvores com 1–3 mudanças de λ (razão 2–4×), com 112 e 250 pontas;
  - comunidades simuladas com γ e μ constantes;
  - n = 14 e n = 40 comunidades.
- **Regimes testados:**
  - (i) os verdadeiros;
  - (ii) os estimados, com DR por ramo (barato), ClaDS (moderado) e BAMM num subconjunto
    de cerca de 50 réplicas (caro).
- **Saída:** proporção de rejeição pelo χ² e pelas partições aleatórias, com IC binomial.
  Meta: 5%.
- Estimativa: 1–2 dias de CPU com 6 núcleos, sem BAMM.

**Resultado S1, testes de razão de verossimilhança (10/11/2026; `data-raw/s1_type1_shifts.R`,
AWS spot, 9000 réplicas em 38 min numa m7i.48xlarge, após duas tentativas interrompidas):**

- Desenho executado: `sim_shift_tree()` com 112 ou 250 pontas, 1 ou 3 mudanças de λ, razão 2
  ou 4, mais controle com razão 1; comunidades com γ = 5 e μ = 1 em todos os ramos; 14 ou 40
  comunidades; 500 réplicas por cenário (18 cenários).
- Regimes verdadeiros, γ e μ livres (`vary = "both"`): rejeição a 5% entre 0,020 e 0,060.
- Regimes verdadeiros, só γ: entre 0,034 e 0,068.
- DR por ramo como covariável de γ: entre 0,036 e 0,066.
- Nenhum IC de Clopper–Pearson fica inteiramente acima de 5%; o teste não é liberal com
  mudanças de taxa na árvore. Com 250 pontas e 3 mudanças, `both` é conservador (0,020 [0,010,
  0,036] e 0,024 [0,012, 0,042] com 14 comunidades).
- Sem convergência: 21 de 9000 réplicas (todas com 3 mudanças); nenhum erro.
- Pendente: regimes estimados (ClaDS, BAMM).

**Resultado S1, teste contra partições aleatórias (10/11/2026; MODE=partition, 50 réplicas por
cenário, 49 partições aleatórias, `vary = "both"`; `data-raw/s1_partition.rds`):**

- **Liberal com 1 mudança de λ.** Rejeição a 5%: 0,12–0,40 nos cenários com 1 mudança
  (média 0,285 com razão 2 ou 4); mesmo nos controles com razão 1 (mudança sem efeito em λ):
  0,12 [0,045, 0,243] com 112 pontas e 0,22 [0,115, 0,360] com 250. Cresce com a razão
  (razão 4: 0,30–0,40).
- Com 3 mudanças: 0,000–0,120 (média 0,055), sem excesso claro.
- O χ² do mesmo ajuste (dentro de `damocles_partition_test()`, 4 inícios) fica em 0,000–0,120
  (média 0,048 com efeito em λ): o problema é a distribuição nula das partições aleatórias,
  não o ajuste.
- Interpretação (inferência, não testada): as partições aleatórias são pareadas pelo número de
  ramos (tolerance 2) e afastadas da verdadeira (max_overlap 0,25), mas o clado com mudança de λ
  tem outra estrutura temporal (ramos curtos, coroa profunda, nasce cedo pelo limiar de
  3 a n/4 linhagens); o LR do clado verdadeiro não é permutável com o de clados aleatórios do
  mesmo tamanho. É justamente o uso pretendido (regimes de BAMM/ClaDS), então o teste, como
  está, não serve para S1–S2 com clados de mudança de λ.
- Caminhos: (a) parear clados aleatórios também por comprimento total de ramos ou idade da
  coroa; (b) trocar o nulo por bootstrap paramétrico (simular comunidades com as taxas
  constantes estimadas e reajustar com a partição verdadeira), permutável por construção;
  (c) reportar só o χ², que está calibrado.
- Custo: o teste de partição levou ~2,5 vezes o extrapolado e sofreu interrupções spot a cada
  10–30 min (23 relançamentos); custo total do S1 estimado em ~US$ 40 (sem acesso ao Cost
  Explorer para confirmar).

**Resultado S1, teste de partição v3 (14/11/2026; depois do PR #76: `random_regimes()` exclui
partições equivalentes à testada e o p conta só partições aleatórias distintas; 99 sorteios,
50 réplicas por cenário, AWS on-demand c7i.8xlarge, run `comdynamite-s1-partition3`;
`data-raw/s1_partition.rds`):**

- **O excesso com 1 mudança sumiu.** Rejeição a 5% nos cenários com 1 mudança: 0,000–0,065
  (antes 0,12–0,40). A interpretação acima (clado com mudança não permutável com clados
  aleatórios) fica superada: boa parte do excesso vinha de partições aleatórias equivalentes à
  verdadeira entrando no nulo.
- **Mas, com 1 mudança, o nulo quase não existe.** Mediana de 3 partições distintas por réplica
  (máximo 47); 79 de 500 réplicas sem nenhuma (sem p), e 325 com menos de 19, que não podem
  rejeitar a 5% (p mínimo 1/(k+1)). A rejeição baixa é, em boa parte, falta de resolução.
- **Com 3 mudanças:** mediana de 93–98 partições distintas; rejeição 0,042–0,191. O cenário 112
  pontas, razão 2, 14 comunidades fica acima de 5% (0,191 [0,091, 0,333]); 112, razão 4, 40
  comunidades, 0,111 [0,037, 0,241].
- **Nas réplicas com ≥ 19 partições distintas** (as que podem rejeitar), a rejeição é 0,083
  [0,037, 0,158] com 1 mudança (n = 96) e 0,081 [0,055, 0,115] com 3 (n = 356): levemente
  liberal.
- O χ² do mesmo ajuste: 0,000–0,120, média 0,047; continua calibrado.
- Conclusão: o teste de partição, como está, não serve para regimes com uma mudança (poucas
  partições comparáveis na árvore) e é levemente liberal quando há partições suficientes.
  Para o artigo, os caminhos (b) bootstrap paramétrico e (c) reportar o χ² ganham peso; o (a)
  não resolve a falta de partições com 1 mudança.
- Execução: 900 réplicas numa c7i.8xlarge on-demand, de 11/11 a 14/11 (~3 dias de máquina);
  nenhum erro.

**Resultado S1, bootstrap paramétrico, rodada reduzida (14/11/2026; PR #99,
`damocles_partition_test(n_rand = 0, n_boot = 49)`; só os 10 cenários com 1 mudança, 50
réplicas cada; spot c7i.48xlarge, run `comdynamite-s1-boot1`; `data-raw/s1_boot.rds`):**

- **Calibrado.** Rejeição a 5% de 0,000 a 0,100 por cenário; agregada, 24/500 = 0,048
  [0,031, 0,071] (IC exato binomial). A distribuição dos p é uniforme (46–55 réplicas por
  decil).
- O χ² das mesmas réplicas: 27/500 = 0,054 [0,036, 0,078]. Os dois concordam.
- Nenhum ajuste falhou (49 de 49 bootstraps válidos em todas as réplicas; zero erros).
- O bootstrap resolve a falta de resolução do teste de partição com 1 mudança: sempre há
  p, com mínimo 1/50.
- Ressalva de interpretação: o bootstrap testa o mesmo nulo do χ² (taxas constantes no
  modelo DAMOCLES), só que com a distribuição do LR simulada em vez da assintótica. Ele não
  isola o efeito do regime de uma má especificação genérica, que era a motivação do teste de
  partição. Aqui o χ² já estava calibrado, então o ganho é de robustez, não de correção.
- Execução: 500 réplicas em ~2 h numa c7i.48xlarge spot, sem interrupções; custo estimado
  em ~US$ 5 (não confirmado no Cost Explorer).
- Os cenários com 3 mudanças foram rodados depois (abaixo). Continua pendente o poder (S2),
  para saber se o bootstrap perde poder em relação ao χ².

**Resultado S1, bootstrap paramétrico, desenho completo (15/11/2026; cenários com 3 mudanças
adicionados, run `comdynamite-s1-boot3`; os 18 cenários estão em `data-raw/s1_boot.rds`):**

| Mudanças | Réplicas | Bootstrap (rejeição a 5%) | χ² (mesmas réplicas) |
|---|---|---|---|
| 1 | 500 | 0,048 [0,031, 0,071] | 0,054 [0,036, 0,078] |
| 3 | 400 | 0,048 [0,029, 0,073] | 0,038 [0,021, 0,061] |

- Com 3 mudanças, a rejeição por cenário vai de 0,00 a 0,08. O cenário em que o teste de
  partição v3 era liberal (112 pontas, razão 2, 14 comunidades: 0,191) dá 0,08 no bootstrap
  e 0,08 no χ², dentro do esperado para n = 50.
- Os p do bootstrap são uniformes nos dois blocos (77–87 réplicas por quintil com 3
  mudanças). Nenhum erro, e 49 de 49 bootstraps válidos em todas as 900 réplicas.
- **Conclusão do S1:** sob diversificação heterogênea e taxas de montagem constantes, o LR
  do DAMOCLES está calibrado tanto pelo χ² quanto pelo bootstrap paramétrico. O teste de
  partição aleatória não serve com 1 mudança e é levemente liberal com 3. Para o artigo:
  χ² como teste principal, bootstrap paramétrico como checagem de robustez, e o teste de
  partição como diagnóstico, não como teste de significância. Isso supõe que o bootstrap não
  perde poder demais, o que o S2 tem de mostrar.
- Execução com 3 mudanças: 400 réplicas, ~14 h de c7i.48xlarge spot em três instâncias
  (duas interrupções spot; uma delas sem relançamento porque o SSO expirou de madrugada).
  Custo estimado em ~US$ 30 (não confirmado no Cost Explorer). Cada réplica com 3 mudanças
  custa ~10 vezes uma com 1 mudança.

### S2. Poder

Desenho revisado em 15/11/2026, depois do S1 (script `data-raw/s2_power.R`):

- Regimes verdadeiros (conhecidos); regimes estimados ficam para o S3. Árvore como no S1
  (`sim_shift_tree()`, uma mudança que dobra λ), sorteada de novo até o clado focal ter a
  fração pedida das pontas (±20%). Taxas de base do S1 (γ = 5, μ = 1); no regime focal, γ, μ
  ou ambos multiplicados pela razão (ambos: a mesma razão, ocupação parecida, mais turnover).
- **Camada A (χ²):** n {5, 10, 14, 20, 40} × pontas {112, 250} × razão {1,1; 1,25; 1,5; 2; 3} ×
  efeito {γ, μ, ambos}, fração 0,25 (150 cenários), mais a fração 0,1 e 0,5 em 112 pontas,
  n {14, 40}, efeito em γ (20 cenários). 200 réplicas; testes p_both (df 2), p_gamma e p_mu (df 1) e a
  razão estimada. O tamanho (razão 1) vem do S1. 50 pontas cortadas (nenhum caso real).
- **Camada B (bootstrap contra χ², mesmas árvores e comunidades):** 112 pontas, n {14, 40},
  as cinco razões, efeito em γ, fração 0,25 (10 cenários); 100 réplicas, 49 simulações.
- Caso de Benício: n = 14 com 112 e 250 pontas, nas duas camadas.

**Piloto (15/11/2026, local, 18 núcleos):** 8 cenários × 6 réplicas do χ² com razão 2 e 2
cenários × 2 réplicas do bootstrap.

- Tempo por réplica: χ² ~5 s com 112 pontas e ~10 s com 250, quase sem efeito de n; bootstrap
  ~40 s em 24 núcleos. Nenhum erro, nenhuma falha de convergência.
- Custo: camada A ~41 h de processo (~2,6 h local com 16 processos; ~15–20 min numa
  c7i.48xlarge); camada B ~7 h local ou <1 h numa 48xlarge. Na AWS spot, as duas juntas ficam
  em ~US$ 3–5 com a subida das instâncias.
- **Com razão 2 e efeito em γ, o poder já é 1 com n = 5** (6/6 réplicas, nos dois tamanhos de
  árvore). Efeito em μ com n = 40: 2–3 de 6. A razão estimada fica perto de 2 (1,84–2,14).
  Com essas taxas, a curva só vai mostrar alguma coisa em razões pequenas e no efeito em μ.
  Por isso as razões 1,1 e 1,25 entraram nas duas camadas (decisão de 15/11/2026): 170
  cenários × 200 réplicas no χ² e 10 × 100 no bootstrap.

**Resultados (16/11/2026, AWS spot; `data-raw/s2_power.rds` e `data-raw/s2_boot.rds`).**
Camada A: 34.000 réplicas, ~177 h de processo (mediana de 12 s por réplica com 112 pontas e
26 s com 250), nenhum erro nem falha de convergência. Camada B: 1.000 réplicas, nenhum erro.
Poder a 5%, fração focal 0,25; linhas n, colunas razão; teste correspondente ao efeito
(γ → p_gamma, μ → p_mu, ambos → p_both).

Efeito em γ:

| n | 112 pontas: 1,1 | 1,25 | 1,5 | 2 | 3 | 250 pontas: 1,1 | 1,25 | 1,5 | 2 | 3 |
|---|---|---|---|---|---|---|---|---|---|---|
| 5 | 0,05 | 0,18 | 0,52 | 0,96 | 1 | 0,11 | 0,36 | 0,86 | 1 | 1 |
| 10 | 0,12 | 0,36 | 0,81 | 1 | 1 | 0,14 | 0,62 | 0,98 | 1 | 1 |
| 14 | 0,14 | 0,48 | 0,95 | 1 | 1 | 0,22 | 0,77 | 1 | 1 | 1 |
| 20 | 0,16 | 0,62 | 1 | 1 | 1 | 0,28 | 0,90 | 1 | 1 | 1 |
| 40 | 0,27 | 0,88 | 1 | 1 | 1 | 0,50 | 0,99 | 1 | 1 | 1 |

Efeito em μ:

| n | 112 pontas: 1,1 | 1,25 | 1,5 | 2 | 3 | 250 pontas: 1,1 | 1,25 | 1,5 | 2 | 3 |
|---|---|---|---|---|---|---|---|---|---|---|
| 5 | 0,04 | 0,04 | 0,05 | 0,10 | 0,18 | 0,04 | 0,05 | 0,05 | 0,12 | 0,30 |
| 10 | 0,04 | 0,06 | 0,10 | 0,12 | 0,44 | 0,04 | 0,06 | 0,10 | 0,20 | 0,52 |
| 14 | 0,06 | 0,06 | 0,09 | 0,16 | 0,50 | 0,04 | 0,06 | 0,06 | 0,28 | 0,60 |
| 20 | 0,06 | 0,05 | 0,10 | 0,26 | 0,62 | 0,06 | 0,06 | 0,09 | 0,36 | 0,74 |
| 40 | 0,06 | 0,09 | 0,16 | 0,39 | 0,87 | 0,04 | 0,10 | 0,20 | 0,49 | 0,92 |

Efeito nos dois (mesma razão), p_both:

| n | 112 pontas: 1,1 | 1,25 | 1,5 | 2 | 3 | 250 pontas: 1,1 | 1,25 | 1,5 | 2 | 3 |
|---|---|---|---|---|---|---|---|---|---|---|
| 5 | 0,04 | 0,11 | 0,34 | 0,76 | 0,99 | 0,04 | 0,22 | 0,74 | 0,97 | 1 |
| 10 | 0,06 | 0,20 | 0,63 | 0,99 | 1 | 0,12 | 0,39 | 0,90 | 1 | 1 |
| 14 | 0,10 | 0,25 | 0,81 | 0,98 | 1 | 0,12 | 0,58 | 0,98 | 1 | 1 |
| 20 | 0,09 | 0,37 | 0,90 | 1 | 1 | 0,15 | 0,72 | 1 | 1 | 1 |
| 40 | 0,15 | 0,72 | 1 | 1 | 1 | 0,32 | 0,94 | 1 | 1 | 1 |

- **γ é bem detectado; μ, não.** Com γ = 5 e μ = 1, a ocupação esperada é γ/(γ + μ) ≈ 0,83, e
  multiplicar μ por 1,5 só a leva a ≈ 0,77. Com n = 14, a razão 1,5 em γ dá 0,95–1, e a razão
  3 em μ dá 0,50–0,60. Árvores maiores ajudam mais que dobrar n entre 5 e 10.
- **O teste livre nos dois (df 2) custa pouco:** p_both fica 0,05–0,10 abaixo do teste de uma
  taxa quando só uma muda.
- **Os testes de uma taxa não atribuem a taxa.** Com efeito só em γ, o teste que libera só μ
  rejeita quase tanto (112 pontas, n = 14: 0,34 com razão 1,25 e 0,87 com 1,5). Com efeito em
  μ, o teste de γ também rejeita (0,43 com razão 3). Um p_gamma ou p_mu significativo indica
  que os regimes diferem, não qual taxa muda; para isso é preciso comparar os modelos γ e μ
  entre si (AIC ou LRT não aninhado), o que ainda não testamos.
- **Fração focal:** com 112 pontas e efeito em γ, a fração 0,1 perde poder (n = 14, razão
  1,25: 0,26 contra 0,48 com 0,25); 0,5 quase não ganha (0,55).
- **Razão estimada sem viés para γ** (mediana/verdadeira entre 0,99 e 1,01). Para μ a mediana
  também fica perto do verdadeiro (0,94–1,03), mas 953 das 10.200 réplicas com efeito em μ
  (428 delas com n = 5) têm μ estimado na fronteira (zero) num dos regimes, com razão < 0,1 ou
  > 10. μ é pouco identificável com poucas comunidades.
- **Bootstrap contra χ² (camada B):** concordam (tabela abaixo; razões 2 e 3 dão 1 nos dois).
  O χ² rejeita sozinho em 26 réplicas e o bootstrap sozinho em 5; com 49 simulações, o p
  mínimo do bootstrap é 0,02, o que explica parte da diferença.

| n | razão | bootstrap | χ² |
|---|---|---|---|
| 14 | 1,1 | 0,15 | 0,13 |
| 40 | 1,1 | 0,14 | 0,21 |
| 14 | 1,25 | 0,35 | 0,41 |
| 40 | 1,25 | 0,68 | 0,76 |
| 14 | 1,5 | 0,86 | 0,88 |
| 40 | 1,5 | 1 | 1 |

- **Caso de Benício (n = 14):** diferenças de 50% em γ entre regimes são detectadas com
  segurança (0,95 com 112 pontas, 1 com 250); de 25%, com 0,48–0,77. Diferenças em μ só
  aparecem se forem de ~3 vezes, e mesmo assim com poder ~0,5–0,6.
**Camada C: qual taxa mudou (16/11/2026, local, 16 processos, 109 min;
`data-raw/s2_compare.rds`, `MODE=compare`).** Os quatro modelos (constante, γ livre, μ livre,
ambos livres) reajustados nas mesmas árvores e comunidades da camada A (mesmas sementes; os
p_gamma conferem exatamente). n {5, 14, 40} × 112 e 250 pontas × as cinco razões × os três
efeitos, fração 0,25, 100 réplicas (9.000), nenhum erro. Os modelos γ e μ têm o mesmo número
de parâmetros, então escolher o de maior verossimilhança é o mesmo que escolher por AIC.

Acerto da taxa entre as réplicas em que o melhor modelo de uma taxa rejeita o constante
(LRT, 5%):

| efeito | pontas | 1,1 | 1,25 | 1,5 | 2 | 3 |
|---|---|---|---|---|---|---|
| γ | 112 | 0,80 | 0,84 | 0,95 | 0,98 | 1 |
| γ | 250 | 0,78 | 0,93 | 1 | 1 | 1 |
| μ | 112 | 0,47 | 0,52 | 0,66 | 0,65 | 0,71 |
| μ | 250 | 0,44 | 0,61 | 0,57 | 0,70 | 0,74 |

- **Efeito em γ detectado é atribuído a γ** em 78–84% das vezes com razão 1,1 e em 95–100%
  com razão ≥ 1,5.
- **Efeito em μ detectado vai para γ em ~30% das vezes,** mesmo com razão 3; com razões
  pequenas o acerto é o do acaso. Sem sinal, a escolha pende para γ (59% das réplicas com
  razão 1,1, somando os três efeitos).
- **Por AIC entre os quatro modelos,** o efeito em μ só é preferido na maioria das réplicas
  com n = 40 e razão 3 (0,65 com 112 pontas, 0,75 com 250); com n = 14 e razão 3, 0,39–0,51.
- **Ambas livres contra o melhor modelo de uma taxa (LRT, df 1):** quando só uma taxa muda,
  rejeita em 0–8% (perto do nominal). Quando as duas mudam, o modelo γ quase sempre vence o
  μ, e a mudança em μ só é detectada com razões grandes: 0,23–0,32 com razão 3 e n = 14,
  0,55–0,69 com razão 3 e n = 40.
- **Comparação com o DDGeoSSE** (Soewongsono & Landis 2026, preprint lido em 16/11/2026):
  lá também o efeito sobre a extinção é o menos preciso, com outro modelo e inferência por
  deep learning; os classificadores deles acertam 73–87% da presença de cada efeito. Os
  desenhos não são comparáveis número a número.

- **Conclusão do S2:** o χ² continua como teste principal (o bootstrap não muda o poder). No
  artigo:
  - reportar o poder por taxa;
  - um teste de uma taxa não atribui o efeito; a atribuição a γ é confiável com efeitos de
    pelo menos 50%, e a atribuição a μ não é (~30% de erro mesmo com razão 3);
  - enquadrar a pergunta empírica como "γ difere entre regimes?" e "a dinâmica local difere
    entre regimes?", e tratar μ como secundário, com intervalos por perfil em vez de
    decisão por teste.
- **Benício como banco de prova (conversa de 16/11/2026):**
  - **Origem das listas (conferido em 16/11/2026, Métodos de Benício et al. 2021):**
    amostragem própria e padronizada, de dezembro de 2014 a fevereiro de 2017, 27 dias por
    área (9 por ano, em três anos), seis sítios reprodutivos por área (duas poças, dois
    riachos, dois transectos em mata), com busca acústica, coleta de girinos e busca visual.
  - **Sobreposição com o ATLANTIC AMPHIBIANS** (sítios a até 15 km de cada área; script em
    scratchpad, a passar para `data-raw/`): todas as 14 áreas têm sítios por perto (o mais
    próximo a 0–14 km; em dez áreas a menos de 2 km). Em quatro
    áreas os registros vêm do mesmo projeto: PEVA (Benício & Da Silva 2017) e PECU, PESB e
    PESV (Da Silva et al. 2017a, Serra do Mar, 2014–2017). Essas quatro ficam fora da
    comparação (circular). Nas outras dez, os estudos do ATLANTIC AMPHIBIANS são
    independentes (1987–2015); seis delas (EECA, EEAS, EEJA, FEENA, EESB, PECB) têm sítio a
    menos de 2 km, e quatro a 4,7–14,1 km (EEIT, EEJI, PETAR, PEJU).
  - Comparar as listas das dez áreas independentes com as dos sítios correspondentes do
    ATLANTIC AMPHIBIANS para estimar a omissão (parametriza a detecção imperfeita).
    Ressalvas: os sítios não são idênticos (distância de corte a decidir), os períodos
    diferem (turnover temporal real se confunde com omissão) e Benício também omite
    espécies, então a comparação dá discordância, não omissão absoluta;
  - **Resultado (16/11/2026, `data-raw/benicio2021_omission.R` e `.rds`).** Fração da
    união das duas listas ausente de cada uma, somando as áreas:

    | corte | áreas | ausente do ATLANTIC AMPHIBIANS | ausente de Benício | Jaccard |
    |---|---|---|---|---|
    | 2 km | 6 | 0,30 (mediana 0,28) | 0,14 (mediana 0,13) | 0,46–0,65 |
    | 15 km | 10 | 0,21 (mediana 0,12) | 0,27 (mediana 0,25) | 0,12–0,81 |

    - Com 2 km (sítios praticamente no mesmo lugar, um ou dois estudos por área), as listas
      do ATLANTIC AMPHIBIANS deixam de fora ~30% das espécies das duas listas juntas, o
      dobro de Benício. O caso extremo é o PECB: 21 espécies só em Benício (muitas de
      riacho e de serapilheira, como *Cycloramphus*, *Ischnocnema*, *Paratelmatobius*),
      contra 8 só no ATLANTIC AMPHIBIANS.
    - Com 15 km, mais estudos entram e a lacuna do ATLANTIC AMPHIBIANS cai para 21%, mas
      as espécies de fora da área inflam o que falta em Benício (27%). O PEJU (um sítio a
      14 km, 13 espécies contra 43) mostra que um sítio distante não serve.
    - Parte da discordância é taxonômica e não omissão: pares como *Scinax similis*/*S.
      fuscomarginatus*, *Physalaemus olfersii*/*P. lateristriga* e *Scinax alter*/*S.
      littoralis*/*S. tymbamirim* aparecem trocados entre as listas da mesma área. Não
      foram resolvidos à mão.
    - **Implicação:** para listas do ATLANTIC AMPHIBIANS, uma taxa de falsa ausência por
      espécie presente de ~0,2–0,3 é um ponto de partida realista para a simulação de
      detecção (variante do S6). Como μ já é o parâmetro mal identificado (S2), esse viés
      pode dominar a estimativa de μ no E6; γ deve ser menos afetado, mas isso precisa ser
      simulado. Filtrar sítios por esforço (2.2) e agregar estudos por localidade são as
      mitigações óbvias a testar;
  - diagnóstico espacial com γ e μ por comunidade contra distância;
  - geografia real das 14 áreas para desenhar o S5.

### S3. Integração do posterior (o ponto metodologicamente aberto)

Hoje `damocles_samples()` resume por "proporção de amostras com p < 0,05" e β mediano. A
diluição já mostrada com o ClaDS (β mediano das amostras 0,22; MAP 1,64; verdadeiro 2)
torna esse resumo inadequado para covariáveis contínuas.

Estimadores a comparar:

- (a) mediana por amostra, o atual;
- (b) taxas médias do posterior por ramo como covariável;
- (c) **verossimilhança integrada sobre configurações**: L̄(θ) = (1/K) Σ_k L(θ | config_k),
  maximizada uma vez. É natural para regimes discretos do BAMM e propaga a incerteza das
  configurações em vez de resumi-la. A ressalva a checar é que as configurações do BAMM
  foram amostradas condicionalmente à árvore, e não aos dados de comunidade;
- (d) correção de erro nas variáveis (SIMEX ou calibração de regressão) para covariáveis
  contínuas, usando a variância posterior por ramo.

Métricas: viés e RMSE de β, cobertura dos intervalos por perfil, erro tipo I e poder.

Entregável provável: um argumento `integrate = TRUE` (ou uma função nova) que implemente
(c), e uma recomendação documentada para (b) ou (d). Isso vira código no pacote, com testes,
num PR próprio.

### S4. Calibração do teste de partições aleatórias

- Indício atual: 0,10 de rejeição no cenário nulo com 40 réplicas e 0,17 num teste menor.
- 500 réplicas do nulo (com e sem diversificação heterogênea) e variação de `tolerance`,
  `max_overlap` e `n_rand`.
- Se o teste for liberal, avaliar:
  - mid-p;
  - restringir o sorteio a clados com idade de coroa parecida, além do tamanho;
  - mais partições por réplica.

### S5. Dependência espacial entre comunidades

- O ajuste compartilhado trata as comunidades como independentes, dada a árvore. Sítios
  próximos compartilham história de colonização, e o LRT tende a ficar anticonservador.
- **Simulação:** comunidades em blocos que compartilham os estados até um tempo de
  separação t_s, e só então evoluem independentemente. t_s = 0 equivale a comunidades
  independentes.
- **Avaliar:**
  - inflação do LRT em função de t_s;
  - bootstrap em blocos de sítios (por região, ou por agrupamento espacial) como correção;
  - partições aleatórias com reamostragem em blocos.
- Com Benício, os blocos naturais são as duas regiões.

#### Diagnóstico em Benício (16/11/2026, `data-raw/s5_benicio_spatial.R` e `.rds`)

- 14 comunidades na árvore de 112 espécies (melhor árvore, altura reescalada para 1), taxas
  compartilhadas, 999 conjuntos de 14 comunidades simuladas independentes com as taxas
  ajustadas. Regiões: interior (latitude ao norte de 23° S, sete áreas) e costa/Serra do Mar
  (sete áreas).

  | estatística | observado | nulo, IC 95% |
  |---|---|---|
  | Jaccard médio entre pares | 0,31 | 0,15–0,18 |
  | Jaccard dentro da região − entre regiões | 0,43 | −0,02–0,02 |
  | Mantel r (Jaccard × distância), todos os pares | −0,62 | −0,21–0,22 |
  | Mantel r dentro das regiões | −0,16 | −0,31–0,32 (p = 0,14) |

- A dependência é regional (as duas regiões quase não compartilham espécies); dentro de cada
  região não há declínio detectável com a distância (só 21 pares por região).
- As taxas diferem entre regiões (LRT 66,9, 2 gl, p ≈ 3 × 10⁻¹⁵). No interior, γ = 3,75 e
  μ = 12,2; na costa, as duas vão ao infinito (presença sem sinal filogenético). O ajuste
  conjunto também vai a taxas enormes (γ ≈ 5000, μ ≈ 12700): com as 14 comunidades e o pool
  de 112, o modelo trata a ocupação como sorteio por espécie.
- Interpretação (minha): filtragem ambiental e biogeografia (interior com Cerrado × floresta
  costeira) somadas ao pool mal especificado (S6), mais do que história de colonização
  compartilhada em escala fina; com esses dados as duas não se separam.

#### Simulação (16/11/2026, `data-raw/s5_spatial.R` e `.rds`)

- Árvores do S2 (112 pontas, regime focal com 0,25 das pontas, γ = 5, μ = 1), 200 réplicas por
  cenário, 4800 no total (55 min, 16 núcleos), sem erros nem falhas de convergência.
  - **História compartilhada:** blocos (n = 14: dois de sete; n = 40: quatro de dez) que
    compartilham estados e resultados de especiação até a idade de separação a_s ∈ {∞, 0,2,
    0,1, 0,05, 0,02}. Correlação dos estados das pontas dentro do bloco ≈ 0,14, 0,25, 0,48 e
    0,7 (0,11 entre blocos).
  - **Regiões com taxas diferentes:** dois grupos de comunidades, μ × 3 ou × 10 no segundo
    (em toda a árvore), mais um teste estratificado (taxas de fundo e efeito de regime livres
    por grupo, verossimilhanças somadas).
- **Tamanho (sem efeito), LRT de γ / μ / ambos a 5%:**

  | a_s | n = 14 | n = 40 |
  |---|---|---|
  | ∞ | 0,06 / 0,05 / 0,05 | 0,03 / 0,05 / 0,05 |
  | 0,2 | 0,09 / 0,08 / 0,08 | 0,07 / 0,07 / 0,06 |
  | 0,1 | 0,12 / 0,12 / 0,12 | 0,22 / 0,28 / 0,27 |
  | 0,05 | 0,34 / 0,29 / 0,39 | 0,32 / 0,41 / 0,46 |
  | 0,02 | 0,47 / 0,36 / 0,58 | 0,44 / 0,44 / 0,60 |

- **História compartilhada infla muito o teste de regime,** e mais com n maior (mais
  comunidades por bloco = mais pseudorreplicação). A estimativa da razão de γ continua sem
  viés (mediana 0,96–1,07 sem efeito; 1,47–1,56 com efeito 1,5). O poder cai pouco
  (n = 14: 0,93 → 0,77 com a_s = 0,02).
- **Taxas diferentes entre regiões não inflam o teste:** 0,03–0,06 com μ × 3 ou × 10, tanto
  no teste conjunto quanto no estratificado (0,03–0,08). A heterogeneidade entre grupos de
  sítios, sozinha, não é o problema.
- **Implicações para o E6:**
  - o perigo é a dependência por história compartilhada entre sítios próximos, não a
    diferença de taxas entre regiões. No ATLANTIC AMPHIBIANS, com centenas de sítios e
    muitos pares a poucos km, a pseudorreplicação é provável; o teste χ² com todas as
    comunidades não é confiável sem correção;
  - correções a avaliar (camada B do S5): (i) agregar sítios em células ou localidades, o
    que também reduz a falsa ausência (S6b); (ii) bootstrap paramétrico com simulação em
    blocos, com a_s estimado da similaridade observada dentro dos blocos; (iii)
    verossimilhança composta com ajuste do LRT por blocos (Chandler & Bate 2007, Biometrika,
    doi:10.1093/biomet/asm015, verificado no Crossref), viável com muitos blocos, inútil com os dois de
    Benício; (iv) bootstrap em blocos espaciais, idem;
  - com Benício (dois blocos), nenhuma correção por reamostragem funciona; o caso só serve
    como ilustração, com o teste estratificado por região como análise principal.

### S6. Pool mal especificado

- Simular com um pool verdadeiro de 250 espécies e ajustar com:
  - (i) o pool verdadeiro;
  - (ii) a união das espécies observadas (como a árvore de 112).
- Medir o viés em γ, μ e β, e o efeito no teste de regime.

### S6b. Detecção imperfeita (feito em 16/11/2026)

- `data-raw/s6_detection.R` e `.rds`: árvores e comunidades como no S2 (112 pontas, regime
  focal com 0,25 das pontas, γ = 5 e μ = 1 no fundo), n = 14 ou 40, efeito nenhum ou razão
  1,5 e 3 em γ ou em μ; 200 réplicas por cenário. Em cada réplica, cada presença vira
  ausência com probabilidade q = 0; 0,1; 0,25; 0,4, sobre a mesma árvore e as mesmas
  comunidades (pareado). Os quatro modelos da camada C do S2 em cada conjunto. Sem erros nem
  falhas de convergência em 8000 ajustes (73 min, 16 núcleos).
- **Viés das taxas absolutas, forte e previsível.** Sem efeito, a mediana de γ cai de 4,98
  (q = 0) para 3,73 (q = 0,25) e 3,00 (q = 0,4), e a de μ sobe de 0,91 para 2,13 e 2,92: a
  falsa ausência vira extinção local. Com q = 0,25, μ sai inflado ~2× e γ ~25% abaixo.
  Taxas absolutas de listas do ATLANTIC AMPHIBIANS não podem ser interpretadas sem correção.
- **O teste de regime continua calibrado.** Sem efeito, a rejeição fica em 0,04–0,09 nos
  três testes em todos os q (sem tendência com q).
- **A razão de γ entre regimes não tem viés:** com efeito em γ, a mediana de γ₁/γ₀ fica em
  1,50–1,54 (razão 1,5) e 2,92–2,98 (razão 3) em todos os q. A razão de μ é atenuada: com
  efeito em μ de razão 3, a mediana cai de 2,90–3,00 (q = 0) para 1,78–1,86 (q = 0,25).
- **Poder cai com q, mais para μ.** Teste de γ com razão 1,5: n = 14 vai de 0,92 (q = 0) a
  0,74 (q = 0,25) e 0,55 (q = 0,4); com n = 40 fica ≥ 0,95. Teste de μ com razão 3: n = 40
  vai de 0,84 a 0,62 (q = 0,25) e 0,51 (q = 0,4); n = 14, de 0,46 a 0,32 e 0,17.
- **Atribuição:** efeitos em γ detectados continuam atribuídos a γ (0,82–1); efeitos em μ
  detectados vão para μ só em 0,42–0,79 das vezes, pior com q alto.
- **Efeito espúrio em μ quando γ muda.** Com efeito forte em γ (razão 3), a falsa ausência
  cria uma diferença aparente de μ entre regimes no modelo com as duas taxas livres (mediana
  de μ₁/μ₀ de 1,6–2,4 com q > 0) e o LRT "duas taxas contra a melhor de uma" rejeita em
  0,17/0,32/0,42 (n = 40, q = 0,1/0,25/0,4), contra 0,04–0,06 sem falsa ausência. Com razão
  1,5 o excesso é pequeno (≤ 0,10). Interpretação (minha, não verificada analiticamente): o
  regime com γ maior tem mais presenças a perder, e o modelo absorve essa perda como
  extinção local maior.
- **Implicações para o E6:**
  - a pergunta central (γ difere entre regimes; "a dinâmica local difere") resiste à
    detecção imperfeita: tamanho calibrado, razão de γ sem viés, poder alto com n = 40;
  - não reportar γ e μ absolutos como taxas biológicas; reportar razões entre regimes;
  - não interpretar "μ também difere" quando γ difere muito: com q ≈ 0,2–0,3 essa conclusão
    é falso positivo em ~1/3 dos casos. Uma diferença em μ só vale com γ semelhante entre
    regimes, ou com correção de detecção;
  - mitigações a testar: agregar estudos por localidade (reduz q), filtrar por esforço, ou
    um modelo de detecção dentro da verossimilhança (extensão de código, não planejada).

### S7. Incerteza de topologia e datação

- Caso empírico: repetir o ajuste principal nos 100 bootstraps datados do Portik e nas
  97 árvores da Feng.
- Simulação: só um teste rápido de que a variação entre árvores não infla o erro tipo I.

### 3.8 Extensões de código que o plano exige

1. Verossimilhança integrada sobre configurações (S3c).
2. Bootstrap em blocos de comunidades (S5), em `damocles_partition_test()` ou numa função
   de reamostragem.
3. Regime que depende da comunidade (costa e interior) nos ramos terminais (2.2a), se a
   opção (a) for escolhida.
4. Funções auxiliares de casamento de nomes e enxerto ficam em `data-raw/`, não no pacote.

## 4. Caso de estudo: ordem das análises

- **E1. Taxonomia e árvore do pool (feito em 2026-11-09):** `data-raw/benicio2021_trees.R`
  constrói as tabelas de 2.1 e 2.2 e uma árvore-mestre por árvore do Portik (melhor e 100
  bootstraps), com todos os enxertos; cada pool é podado dela. Idades da raiz: 115 Ma para
  112 e para o pool regional (bootstraps de 106 a 154) e 166 Ma para 248 (160 a 195). A
  distância RF normalizada à árvore de Benício é 0,19 para 112 e 0,20 para 248.
  - **Piloto de tempo** (melhor árvore, 14 comunidades, um núcleo): `damocles_fit(shared =
    TRUE)` leva 0,8 s com 112 espécies e 1,6 s com 248; com DR como covariável de γ, 2,4 s e
    3,3 s; um ajuste por comunidade, 10 s e 14 s. Com 6 núcleos, 1.000 réplicas de S1 com
    três modelos cada levam algo como 15–30 min por cenário, e os 100 bootstraps, poucos
    minutos. Com o pool regional (222), o ajuste constante leva 1,4 s. O custo não limita o
    plano com estes tamanhos; pools de milhares de espécies
    (E6) ainda precisam de piloto próprio.
- **E2. Regimes e taxas:**
  - BAMM do Vasquez-Restrepo: o zip tem o event data de cada uma das 100 corridas de
    anfíbios (conferido em 2026-11-09; baixado só o que é de anfíbios);
  - ClaDS próprio sobre o Portik podado (RPANDA, em `data-raw/`);
  - DR com `dr_rates()`;
  - GeoSSE de Benício et al. (2021) como referência da diferença costa contra interior.
- **E3. Ajuste base:** `damocles_fit(shared = TRUE)` com γ e μ constantes, pools de 248 e
  112, e `sympatry` livre e fixa em 0.
- **E4. Regimes e covariáveis:**
  - LRT e intervalos por perfil;
  - `damocles_samples()` sobre as configurações do posterior, com o estimador vencedor
    de S3;
  - `damocles_partition_test()` com a melhor configuração.
- **E5. Robustez:**
  - bootstrap em blocos (S5);
  - 100 árvores do Portik e 97 da Feng (S7);
  - pool de 112 contra 248;
  - com e sem as espécies enxertadas.
- **E6. Escala da Mata Atlântica (parte do artigo único):** o ATLANTIC AMPHIBIANS filtrado por esforço e por qualidade
  da coordenada, com ajustes por clado (famílias ou gêneros com pelo menos cerca de 20
  espécies) e regimes dentro de cada clado. Comparar com os modelos "global" e "por clado"
  de Reijenga et al. (2026).

## 5. Cronograma (estimativa)

| Etapa | Duração |
|---|---|
| E1 (taxonomia e árvore) + piloto de tempo de computação | 1 semana |
| S1, S2, S4, S6 (em paralelo, em background) | 3–5 semanas |
| S3 + código da verossimilhança integrada + S5 e bootstrap em blocos | 2–3 semanas |
| E2–E5 (caso de estudo) | 2–3 semanas |
| E6 (escala da Mata Atlântica, com piloto de tempo para pools grandes) | 3–5 semanas |
| Redação do artigo único | 6–8 semanas |

## 6. Decisões

Tomadas em 2026-11-09:

1. Artigo único para a AmNat.
2. O pool de 248 é amostragem filogenética (sequências no GenBank), não pool regional.
3. Baixar a parte de anfíbios do Figshare, depois de conferir o conteúdo.
4. Contatar Etienne e Pigot só quando houver resultado concreto.

5. Pool regional: Mata Atlântica da Serra do Mar e do interior, com o Cerrado (2.2).
6. *Ischnocnema guentheri* C1, C3 e C4 como *Ischnocnema* sp1–sp3 (2.1).

Pendentes: nenhuma no momento.

## 7. Riscos

- **Poder com n = 14 pode ser baixo** (S2). Com o artigo único, a inferência já se apoia
  também no ATLANTIC AMPHIBIANS (E6), e Benício entra como o caso de amostragem padronizada.
- **Super-árvores com espécies imputadas** (Vasquez-Restrepo) distorcem taxas perto das
  pontas. Mitigação: comparar com o Portik, que é molecular.
- **Confusão entre barreira biogeográfica e dinâmica local** quando costa e interior
  compartilham o pool (2.2).
- **Esforço heterogêneo no ATLANTIC AMPHIBIANS:** ausências falsas inflam μ. A seleção de
  estudos do próprio conjunto (só inventários de comunidades, com composição e coordenadas)
  reduz o problema (2.2), mas não fixa esforço mínimo. Pode exigir
  um modelo de detecção. Com E6 dentro do artigo único, esse risco fica central: filtrar
  sítios por esforço e simular o efeito de detecção imperfeita sobre γ e μ (variante de S6).

## Referências (verificadas no Crossref, OpenAlex ou servidor de preprints em 2026-11-09)

- Benício, R., Provete, D. B., Lyra, M., Heino, J., Haddad, C. F. B., Rossa-Feres, D. C., &
  da Silva, F. R. (2021). Differential speciation rates, colonization time and niche
  conservatism affect community assembly across adjacent biogeographical regions.
  *Journal of Biogeography*, 48, 2211–2225. doi:10.1111/jbi.14145
- Pigot, A. L., & Etienne, R. S. (2015). A new dynamic null model for phylogenetic
  community structure. *Ecology Letters*, 18, 153–163. doi:10.1111/ele.12395
- Portik, D. M., Streicher, J. W., Blackburn, D. C., Moen, D. S., et al. (2023). Redefining
  possible: combining phylogenomic and supersparse data in frogs. *Molecular Biology and
  Evolution*, 40, msad109. doi:10.1093/molbev/msad109
- Reijenga, B. R., Etienne, R. S., Murrell, D. J., & Pigot, A. L. (2026). Speciation
  history shapes patterns of assemblage species richness in birds. *Ecology Letters*, 29,
  e70432. doi:10.1111/ele.70432
- Vancine, M. H., Duarte, K. S., de Souza, Y. S., et al. (2018). ATLANTIC AMPHIBIANS: a
  data set of amphibian communities from the Atlantic Forests of South America. *Ecology*,
  99, 1692. doi:10.1002/ecy.2392
- Vasquez-Restrepo, J. D. (2024). A dataset containing speciation rates, uncertainty, and
  presence-absence matrices for more than 34,000 vertebrate species. bioRxiv (preprint,
  sem revisão por pares). doi:10.1101/2024.04.09.588748
- Para o restante do contexto metodológico (STRAPP, ES-sim, Title & Rabosky 2019, ClaDS,
  BAMM), ver a documentação de `damocles_samples()` e `dr_rates()`.
