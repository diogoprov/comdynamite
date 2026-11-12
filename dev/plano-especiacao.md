# Plano: especiação acoplada à metacomunidade

Objetivo: fechar o limite registrado em `dev/notas-artigo-taxas.md` (ponte Leibold–Januario).
Hoje o comdynamite simula a metacomunidade sobre um pool fixo de espécies e uma árvore dada;
a especiação não depende da dinâmica. Com populações locais que divergem e viram espécies, a
árvore passa a emergir da metacomunidade, e dá para perguntar com verdade conhecida se a
persistência metapopulacional prediz a especiação (Januario et al. 2025) e se a montagem
deixa rastro nas taxas de evolução de traits entre espécies (Leibold et al. 2019, 2022).

Referências conferidas no Crossref ao fim deste arquivo. O que não é afirmação de uma fonte é
proposta nossa.

## Situação atual do código

- `mc_dynamics()` (em `R/dynamics.R`) integra abundâncias de k espécies em m sítios com
  crescimento dependente do ambiente (`r`), competição (`alpha`), dispersão, imigração do
  pool e extirpação. k é fixo: as matrizes têm uma linha por espécie desde o início.
- Os traits são por espécie, não por população; não há evolução durante a dinâmica.
- `trait_mk()`, `history_rates()` e `trait_rates()` simulam traits e taxas sobre uma árvore
  dada. `food_web_dynamics()` tem evolução de traits médios por gradiente, mas sem espaço.

## Modelo proposto

Segue a família de modelos mecanísticos eco-evolutivos em que especiação, dispersão e
seleção ambiental geram a diversidade regional, revisada por Pontarp et al. (2019) e
implementada de forma geral no gen3sis (Hagen et al. 2021). O que o comdynamite acrescenta é a
ligação com os métodos de metacomunidade e de taxas por ramo do pacote, com a verdade
conhecida em cada etapa.

### 1. Populações com trait local

- Unidade: população = espécie × sítio, com abundância (já existe) e trait médio local z_is
  (novo).
- Seleção: z_is anda pelo gradiente de ln W no sítio, como em `food_web_dynamics()`
  (dz/dt = V ∂ln W/∂z, McPeek 2017, eq. 3.7), com W dado pelo ajuste ao ambiente local e,
  opcionalmente, pela competição dependente do trait.
- Fluxo gênico: depois da dispersão, z_is vira a média ponderada pelas abundâncias de
  residentes e imigrantes. Com dispersão alta isso afoga a adaptação local (gene swamping,
  Leibold et al. 2022); com dispersão baixa permite os efeitos de prioridade mediados pela
  evolução de Leibold et al. (2019).
- Variância genética V fixa por espécie na primeira versão; genética explícita (locos, como no
  modelo baseado em indivíduos de Leibold et al. 2019) fica como extensão.

#### Sugestões vindas do ECY25-0086 (2026-11-11; 1, 3, 4, 6 e 8 aprovadas para a Fase 1 em 2026-11-12)

O planejamento de um ABC para separar deriva, seleção e plasticidade na morfologia de
*Dendropsophus nanus* (21 lagoas de Campo Grande) concluiu que o comdynamite atual não serve:
`local_adaptation`, `founder_sd` e `monopolization` são desvios fenomenológicos aplicados em
`sample_individuals()`, sem h², G nem gerações, e plasticidade e adaptação são o mesmo
parâmetro. Plano completo do projeto, fora deste repositório:
`[2023] Cap 1 Marcos/docs/ibm_simulation_plan.md`, seções 2 e 4.3. O que tornaria a seção 1
útil para esse tipo de pergunta:

1. Deriva explícita: z_is recebe N(0, V/Ne_is) por geração, com Ne_is ligado a N_is por um
   parâmetro Ne/N. Sem isso não há nulo neutro dentro do modelo, só seleção e fluxo gênico.
2. Versão multivariada: G por espécie (aproveitando as covariâncias de `traits_mv`), resposta
   de Lande Δz̄ = G (P + Ω)⁻¹ (θ − z̄) com aptidão gaussiana e P = G + E, com E não
   proporcional a G. Permite avaliar o teste B ∝ W (Ackermann & Cheverud 2004) quando W ≠ G,
   cujo erro tipo I Prôa et al. (2013) mediram por simulação. É o item mais caro; pode vir
   depois.
3. Plasticidade separada da genética: fenótipo = média genética + b·E_natal + ruído. A parte
   plástica não é herdada e os migrantes carregam o fenótipo da origem (desenvolvimento
   larval). A saída separa a média genética e a plástica por população, a única forma de ter
   verdade conhecida para seleção contra plasticidade.
4. Fundador mecanístico: população nova começa com média sorteada de N(z̄_fonte, V/k), com k
   fundadores ponderados pela dispersão. É a versão de processo do `founder_sd`, depois
   erodida pelo fluxo gênico.
5. Ótimo que segue um ambiente variável no tempo por sítio (por exemplo, séries anuais de
   cobertura urbana), com registro do tempo desde a mudança. A defasagem da resposta (seleção)
   contra a resposta imediata (plasticidade) é a principal assinatura que separa os dois.
6. Parâmetro explícito de passos por geração (1–2 anos em anuros), para converter a mudança
   em haldanes.
7. Amostragem alocrônica: `sample_individuals()` em vários `record_steps`, para simular séries
   de museu além do retrato espacial.
8. Testes com verdade analítica: (a) população isolada e neutra, Var(Δz̄) = h²t/Ne em
   unidades de desvio padrão (Lande 1976); (b) equilíbrio migração–seleção, divergência da
   ordem de s/(s + m); (c) B ∝ W com inclinação esperada 1 quando G ∝ P.

Ressalva do próprio projeto: com um único retrato espacial, seleção e plasticidade para o
mesmo ótimo ficaram quase inidentificáveis (urbanização atual e acumulada com r = 0,985). Para
esse caso, o valor do bloco está em simular também dados temporais (itens 5 e 7). Avaliação
nossa: os itens 1, 3, 4 e 8 cabem na Fase 1; 5 e 7 tornam o bloco útil para o ABC; 2 e 6 são
extensões.

#### Decisões da Fase 1 (Diogo, 2026-11-12)

- A evolução entra como argumento `evolution` de `mc_dynamics()` e `simulate_metacommunity()`,
  com as opções em `evolution_control()`; `NULL` por padrão, e com `genetic_var = 0` as
  abundâncias são as mesmas de antes para a mesma semente. A função nova da Fase 2 reaproveita
  esse motor.
- Entram as sugestões 1 (deriva, `ne_ratio`), 3 (plasticidade como norma de reação linear em
  torno do ambiente médio, com a parte plástica separada na saída), 4 (amostragem de
  fundadores), 6 (`steps_per_generation`) e 8 (testes analíticos: deriva contra Lande 1976 e
  equilíbrio migração–seleção de duas populações). A 2 (G multivariada), a 5 (ótimo variável
  por sítio, já possível com `env_time`, sem o registro do tempo desde a mudança) e a 7
  (amostragem alocrônica) ficam para depois.
- O crescimento usa o fenótipo médio e ignora a carga da variância em torno dele, para que
  `genetic_var = 0` reproduza o modelo sem evolução.

### 2. Divergência e especiação

#### Decisões da Fase 2 (Diogo, 2026-11-12)

- Nome: `simulate_evolving_metacommunity()`. Modos na primeira versão: `"geographic"`,
  `"ecological"`, `"point"`, `"fission"` e a camada de conclusão prolongada.
- Dois PRs: 2a (função, registro filogenético, extinção global, `"geographic"` e
  `"ecological"`, testes de consistência e validação 3) e 2b (conclusão prolongada, `"point"`,
  `"fission"` e as validações 1 e 2).
- `mc_dynamics()` refatorada em preparo, passo e fechamento sobre um mesmo estado, que a função
  nova usa; sem especiação, ela reproduz `simulate_metacommunity()` com a mesma semente.
- Divergência geográfica com duas regras: clusters de fluxo (padrão; diverge entre clusters,
  zera dentro) e gradual (cresce 1 por geração e cai `homogenization` por migrante).
- Competição que decai com a distância patrística, com piso `inter_min` (sem ele, espécies
  distantes não competem e a riqueza explode).
- Estado da 2a (2026-11-12): implementada e testada; validação 3 sem apoio (ρ ≈ 0 entre
  conectividade e endemismo; `data-raw/evolving_peripheral.R`).

#### Revisão da 2a antes da 2b (Diogo, 2026-11-12)

1. Modo ecológico: persistência mais fluxo reduzido. Um par de populações conta as gerações
   seguidas com |Δz̄| ≥ limiar e taxa de migração < `migration_max` (padrão 0,001); os grupos
   se separam quando a contagem chega a `min_generations` (padrão 20). Sem isso, deriva e
   amostragem de fundadores em populações pequenas especiavam num ambiente uniforme.
2. Regra de clusters: ligação pela taxa de migração (migrantes esperados na geração sobre o
   tamanho da população receptora, o maior dos dois sentidos) ≥ `migration_min` (padrão
   0,005, calibrado em `data-raw/evolving_migration_min.R`: nada especia com 0,002 ou menos, e com 0,02 ou mais quase nenhuma ligação se mantém). Dentro de um cluster a divergência fica congelada e só zera depois de
   `reset_generations` (padrão 10) gerações seguidas no mesmo cluster.
3. Brotamento: o maior grupo (mais sítios, depois mais indivíduos) mantém a identidade e a
   idade da mãe; os outros viram espécies novas. Vários de uma vez saem em bifurcações
   sucessivas, a mais divergente primeiro, com ramos internos de comprimento zero.
4. Competição pelo nicho: α_ij(s) também cai com a sobreposição dos nichos gaussianos dos
   fenótipos no sítio, exp(−(p_is − p_js)²/(2(σ_i² + σ_j²))) (`niche_overlap = TRUE`).
5. `inter_min = 0,25` mantido; sem fusão em contato secundário.

#### Decisões da Fase 2b (Diogo, 2026-11-12)

- `"fission"`: as duas partições, por argumento (`partition = "spatial"`, padrão: as populações
  mais próximas, pelo kernel, de um sítio sorteado, em número sorteado; `"random"`: cada
  população com probabilidade 1/2). Taxa por espécie por geração.
- `"point"`: probabilidade por indivíduo por geração; a espécie nova começa com 1 indivíduo no
  sítio e herda a média genética local, com salto gaussiano opcional (`jump_var`, padrão 0),
  que vale para todos os modos.
- Conclusão prolongada: as incipientes são linhagens independentes na dinâmica (como no PBD);
  a conclusão é só registro. No fim, cada linhagem viva entra na espécie boa mais próxima
  (ela mesma ou um ancestral), representada na árvore pela própria espécie boa ou, se extinta,
  pela incipiente mais antiga.
- Validações 1 (point contra Ewens, θ = 2Jν) e 2 (γ e ajuste do PBD, GPL, só em data-raw).


Hubert et al. (2015) organizam os modelos de especiação em metacomunidades em três quadros:
especiação geográfica (limitação de dispersão e estocasticidade), especiação ecológica
baseada em nicho (adaptação com fluxo gênico em paisagens heterogêneas) e especiação
ecológica mediada pela comunidade (trade-offs de história de vida em paisagens homogêneas).
Também destacam que a especiação tem duas etapas: espécies incipientes surgem continuamente
(a taxa real) e só viram espécies boas se sobreviverem a um período de transição (a taxa
aparente, a que se observa), ideia formalizada na especiação prolongada (Rosindell et al.
2010; Etienne & Rosindell 2012). Os modos abaixo seguem esse quadro; a regra de especiação
condiciona os resultados e precisa ser variada nos experimentos.

**Modos de iniciação** (o que cria uma espécie incipiente):

- `"point"`: mutação pontual; cada nascimento vira espécie nova com probabilidade ν. É o
  modo do MESS (Overcast et al. 2021) e da teoria neutra (Hubbell 2001, como citado por
  Hubert et al. 2015). Serve de limite neutro.
- `"fission"`: fissão aleatória; uma espécie é dividida em duas ao acaso (Hubbell 2001, como
  citado por Hubert et al. 2015; também opção do MESS). Aqui a divisão é de populações, não
  de indivíduos, o que dá às filhas distribuições geográficas próprias.
- `"geographic"` (modo principal): uma matriz de divergência entre populações coespecíficas
  cresce sem fluxo gênico e é reduzida quando há fluxo; o grupo de populações que passa do
  limiar se separa. É a abordagem do gen3sis (Hagen et al. 2021). Variante a considerar: a
  acumulação de incompatibilidades por mutação e deriva de Desjardins-Proulx & Gravel (2012a,
  b), cuja previsão (comunidades periféricas e pouco conectadas acumulam mais espécies por
  especiação) vira um teste qualitativo.
- `"ecological"`: separação quando a diferença de trait entre grupos de populações passa de
  um limiar, mesmo com fluxo gênico (especiação ecológica baseada em nicho).
- Mediada pela comunidade (trade-offs competição–colonização etc.): fica para depois, porque
  exige traits de história de vida que o simulador ainda não tem.

**Conclusão da especiação** (camada sobre qualquer modo, em vez de um modo à parte):

- Cada iniciação cria uma espécie incipiente, que vira espécie boa depois de um tempo τ
  (fixo ou exponencial com taxa λ2, como no modelo PBD) se não se extinguir nem se fundir de
  novo com a mãe (no modo `"geographic"`, fluxo gênico que zera a divergência).
- `completion = 0` dá especiação instantânea. No modo `"geographic"` o limiar de divergência
  já é uma forma mecanística de especiação prolongada; τ adiciona o período de transição de
  Hubert et al. (2015).
- Registrar os tempos de iniciação e de conclusão. A árvore reconstruída sai com os nós nos
  tempos de iniciação, como no PBD (Etienne & Rosindell 2012), e, opcionalmente, nos de
  conclusão. As taxas real (iniciação) e aparente (conclusão) ficam disponíveis, o que importa
  para o teste de Januario, que usa taxas nas pontas.

**Trait das espécies novas:**

- Mecanístico (padrão): as filhas herdam as populações e os traits locais que já tinham.
- Pontuado (opção): salto gaussiano no momento da especiação, como no MESS, que usa variância
  σ²_M/(λ + λε), a variância esperada entre mãe e filha no BM da metacomunidade (Overcast et
  al. 2021). Útil nos modos `"point"` e `"fission"`, em que a divisão não vem de divergência.

**Primeira versão (proposta):** `"point"` e `"fission"` (para validar contra a teoria neutra e
o PBD), `"geographic"` com a camada de conclusão e `"ecological"`.

### 3. Registro filogenético e da história

- Cada especiação registra tempo, espécie mãe, espécies filhas e as populações de cada uma.
  No fim, montar um `phylo` completo (com extintas) e o reconstruído (só sobreviventes, via
  `ape::drop.fossil()` ou equivalente).
- Extinção global quando a abundância regional de uma espécie zera.
- A distribuição de cada linhagem ao longo do tempo (sítios ou tipos de hábitat ocupados) vira
  um `simmap`, pronto para `history_rates()` e `as_herodotools()`.
- Séries por espécie: ocupação, abundância regional e as duas métricas de persistência de
  Januario et al. (2025), ocupação de equilíbrio e resiliência (conferir a definição exata no
  artigo antes de implementar).
- Traits por população e médias por espécie, para os estimadores de taxa por ramo e para as
  estatísticas T (`cati::Tstats()`).

### 4. Interface (decidida em 2026-11-08: função nova)

- Opção A: função nova, por exemplo `simulate_diversification(landscape, ...)`, que reaproveita
  as partes internas de `mc_dynamics()` (dispersão, imigração, competição) e devolve uma
  classe própria com `tree`, `tree_full`, `history` (simmap), `abundance`, `pop_traits`,
  `events` e `persistence`.
- Opção B: argumento `speciation = ` em `simulate_metacommunity()`. Mais simples para o
  usuário, mas complica uma função que hoje assume k fixo.
- Decisão (Diogo, 2026-11-08): opção A, função nova, com um `as_metacommunity()` que converta
  a saída para os métodos que já existem. O nome final fica para a rodada de UX.

### 5. Implementação

- R base, sem dependência nova. O número de espécies cresce: guardar abundâncias e traits em
  matrizes que crescem por blocos (dobrar a capacidade) e um teto `max_species` com erro
  claro.
- Medir o custo cedo (como em `simulation_cost()`), porque a matriz de divergência é
  quadrática no número de populações por espécie. Se ficar lento demais, discutir Rcpp antes
  de adotar, porque é dependência nova.

## Validação

1. Limite neutro (`"point"`, sem seleção): o número de espécies e as abundâncias devem bater
   com as expectativas da teoria neutra para o θ correspondente. Antes, levantar quais
   expressões analíticas usar (Rosindell et al. 2010 dá resultados para o caso prolongado).
2. Modo prolongado contra o PBD: as árvores reconstruídas devem ter a desaceleração aparente
   descrita por Etienne & Rosindell (2012). Ajustar o PBD (pacote `PBD`, conferir a licença) em
   `data-raw` e ver se recupera as taxas.
3. Desjardins-Proulx & Gravel (2012a): no modo `"geographic"` em grafos de sítios, as
   comunidades periféricas e pouco conectadas devem acumular mais espécies por especiação.
4. gen3sis (CRAN, GPL; só em `data-raw`): mesma paisagem e parâmetros equivalentes no modo
   alopátrico; comparar riqueza, tamanho de clado e forma da árvore.
5. Leibold et al. (2019): sem especiação e com duas espécies, reproduzir a neutralidade
   regional qualitativamente (convergência das distribuições de trait e das ocupações) nos
   eixos de extinção de manchas e dispersão da Fig. 5.
6. MESS (Overcast et al. 2021; Python, sem licença declarada no GitHub, só executado): no
   modo `"point"` com um só sítio e imigração de um pool fixo, comparar riqueza e
   abundâncias locais.
7. Consistência interna: árvore ultramétrica no reconstruído, número de eventos igual ao de
   nós, simmap coerente com as ocupações, testes de que sem especiação o resultado é o de
   `mc_dynamics()` com a mesma semente.

## Experimentos que isso abre

- Teste de Januario com verdade conhecida: persistência por espécie contra a taxa de
  especiação nas pontas da árvore reconstruída (por exemplo a estatística DR de Jetz et al.
  2012), variando dispersão e extinção de manchas. Pergunta: o desacoplamento aparece mesmo
  quando a persistência causa a especiação no modelo?
- Neutralidade regional e taxas de trait: taxas por ramo verdadeiras e estimadas com médias de
  espécies, contra a divergência entre populações e as estatísticas T.
- Escalas dos dados: Overcast et al. (2021) mostram que abundâncias, diversidade genética e
  traits/filogenia guardam informação de escalas de tempo diferentes. As saídas por população
  e por espécie permitem perguntar em que eixo o sinal da montagem sobrevive.
- Ligação com os cenários 1, 2, 5 e 8 das notas, agora com a árvore gerada pela própria
  metacomunidade.

## Fases

1. Traits por população com seleção e fluxo gênico em `mc_dynamics()`, sem especiação;
   validação 5.
2. Divergência, os modos de especiação, a camada de conclusão e o registro filogenético;
   validações 1, 2, 3 e 7.
3. Saídas (simmap, persistência, traits por população) e integração com `history_rates()`,
   `trait_rates()` e as estatísticas T.
4. Comparação com o gen3sis e o MESS (validações 4 e 6) e medida de custo.
5. Experimentos e o artigo.

## Riscos e decisões em aberto

- A regra de especiação decide muito do resultado; por isso os vários modos e a obrigação de
  reportar a sensibilidade.
- Escalas: a dinâmica roda em passos ecológicos e a especiação precisa de muitos passos.
  Januario et al. (2025) atribuem o desacoplamento justamente à diferença entre décadas e
  milhões de anos; o modelo comprime essas escalas e o artigo tem de dizer isso.
- Custo computacional com muitas espécies e sítios.
- Decisões para o Diogo: nome da função nova; confirmar os modos da primeira versão; número da versão (proposta: bloco próprio, número a definir, depois da v0.7.0).

## Referências

- Ackermann, R. R., & Cheverud, J. M. (2004). Detecting genetic drift versus selection in human evolution. *PNAS*, 101, 17946–17951. https://doi.org/10.1073/pnas.0405919102
- Desjardins-Proulx, P., & Gravel, D. (2012a). *The American Naturalist*, 179, 137–144.
  https://doi.org/10.1086/663196
- Desjardins-Proulx, P., & Gravel, D. (2012b). *Ecology and Evolution*, 2, 1781–1790.
  https://doi.org/10.1002/ece3.292
- Etienne, R. S., & Rosindell, J. (2012). *Systematic Biology*, 61, 204–213.
  https://doi.org/10.1093/sysbio/syr091
- Hagen, O., Flück, B., Fopp, F., Cabral, J. S., Hartig, F., Pontarp, M., Rangel, T. F., &
  Pellissier, L. (2021). *PLOS Biology*, 19, e3001340.
  https://doi.org/10.1371/journal.pbio.3001340
- Januario, M., Pinsky, M. L., & Rabosky, D. L. (2025). *Ecology Letters*, 28, e70021.
  https://doi.org/10.1111/ele.70021
- Hubert, N., Calcagno, V., Etienne, R. S., & Mouquet, N. (2015). *Ecology Letters*, 18,
  864–881. https://doi.org/10.1111/ele.12458
- Jetz, W., Thomas, G. H., Joy, J. B., Hartmann, K., & Mooers, A. O. (2012). *Nature*, 491,
  444–448. https://doi.org/10.1038/nature11631
- Lande, R. (1976). Natural selection and random genetic drift in phenotypic evolution. *Evolution*, 30, 314–334. https://doi.org/10.1111/j.1558-5646.1976.tb00911.x
- Leibold, M. A., Urban, M. C., De Meester, L., Klausmeier, C. A., & Vanoverbeke, J. (2019).
  *PNAS*, 116, 2612–2617. https://doi.org/10.1073/pnas.1808615116
- Leibold, M. A., Govaert, L., Loeuille, N., De Meester, L., & Urban, M. C. (2022). *Annual
  Review of Ecology, Evolution, and Systematics*, 53, 299–326.
  https://doi.org/10.1146/annurev-ecolsys-102220-024934
- McPeek, M. A. (2017). *Evolutionary Community Ecology*. Princeton University Press.
  https://doi.org/10.23943/princeton/9780691088778.001.0001
- Overcast, I., Ruffley, M., Rosindell, J., Harmon, L., et al. (2021). *Molecular Ecology
  Resources*, 21, 2782–2800. https://doi.org/10.1111/1755-0998.13514
- Pontarp, M., et al. (2019). *Trends in Ecology & Evolution*, 34, 211–223.
  https://doi.org/10.1016/j.tree.2018.11.009
- Prôa, M., O'Higgins, P., & Monteiro, L. R. (2013). Type I error rates for testing genetic drift with phenotypic covariance matrices: a simulation study. *Evolution*, 67, 185–195. https://doi.org/10.1111/j.1558-5646.2012.01746.x
- Rosindell, J., Cornell, S. J., Hubbell, S. P., & Etienne, R. S. (2010). *Ecology Letters*,
  13, 716–727. https://doi.org/10.1111/j.1461-0248.2010.01463.x
