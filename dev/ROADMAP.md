# comdynamite: estado e próximos passos

Atualizado em 2026-11-07. Versão atual: 0.6.2, com `0.6.2.9000` em desenvolvimento (DOI de
conceito `10.5281/zenodo.23993533`).
Tags: `v0.1.0`, `v0.2.0`, `v0.5.0`, `v0.5.1`, com o que faltava da v0.5.0, `v0.5.2`, `v0.6.0`,
`v0.6.1` e `v0.6.2`. O marco v0.3.0 fechou em 2026-11-07 e saiu na 0.6.2, sem tag `v0.3.0`. A
`v0.4.0` começou com `mc_data()` e traits externos e tem 3 itens abertos. Todas as pendências das etapas
anteriores estão fechadas.
As versões abaixo são uma sugestão: cada marco tem um critério de conclusão, e a
ordem pode mudar. Tag no formato `vX.Y.Z`, com release no GitHub.

## Visão

O comdynamite deixa de ser só um simulador e passa a ser um hub de métodos que ligam
traits, filogenias, biogeografia e metacomunidades. A simulação continua no centro:
cada método do pacote pode ser testado com dados em que a verdade é conhecida.
Métodos publicados que nunca chegaram ao R, como os de Leibold et al. (2010) e
Peres-Neto et al. (2012), são pouco usados por isso e têm prioridade.

Para isso funcionar, o pacote precisa de uma estrutura de dados comum, dependências
leves (pacotes grandes em `Suggests`) e convenção de nomes (`simulate_*` para gerar,
e um prefixo a definir para analisar). O título e o `Description` mudaram na 0.6.0.

## Feito

| Item | O que entrou |
|---|---|
| Modelo | Thompson et al. (2020) com traits e competição filogenéticos; paisagens implícitas e em grade; permutações do toro; ambiente variável no tempo |
| Montagem por regras | `assemble_community()`: similaridade limitante, filtragem e neutra |
| Traits | `traits_ou()` (um ou vários ótimos) e `paint_regimes()`, no lugar de `traits_bm()` |
| Identidade | hex sticker (opção F1) e `data-raw/hex_sticker.R` |
| Site | pkgdown em https://provetelab.org/comdynamite/, com build automático no `main` |
| Simulador estático | `simulate_static()`, na versão de desenvolvimento |
| Vignettes | `comdynamite`, `phylogenetic-structure` e `evaluating-methods` |
| Documentação | `\value` e exemplos em todas as funções exportadas; `Description` atualizada |
| Desde a 0.6.1 | `pcps_test()`, toro retangular, `summary_stats()`, `run_scenarios()`, `equilibrium_check()`, saída longa, paralelismo opcional, nulo restrito por conjunto em `trait_partition()`, piloto de identificabilidade |

## v0.1.0: primeira versão citável (concluída em 2026-11-05)

CI em macOS, Windows e Linux; 100% de cobertura de linhas; `inst/CITATION` e
`CITATION.cff`; release no GitHub com DOI do Zenodo.

- [x] Submissão ao r-universe: `comdynamite` incluído em `packages.json`, no repositório
      `diogoprov/diogoprov.r-universe.dev`.
- [x] Build confirmado em https://diogoprov.r-universe.dev/comdynamite (0.6.1, OK em
      Linux, macOS, Windows e wasm); badge e instalação pelo r-universe no README.

## v0.2.0: validação e geradores estáticos

Objetivo: mostrar que a implementação reproduz o modelo publicado e ter uma linha de
base sem dinâmica.

- [x] `simulate_static()`: resposta gaussiana de Dray & Legendre (2008), usada também
      por Peres-Neto et al. (2012). Nessas fontes, `h` vem de U(0,5; 1), a abundância é
      o valor esperado e a incidência é binomial com probabilidade `y`, opções de
      `sampling`. As contagens de Poisson e `h` em U(0, 30), padrão no pcps-sim, vêm de
      Peres-Neto & Kembel (2015, Appendix A) e Duarte et al. (2016, Appendix S2). Duarte
      et al. atribuem o Poisson a 2012, mas ele não está no texto principal desse artigo.
- [x] Dois gradientes em `simulate_static()`: produto das respostas ao ambiente e ao
      espaço, como em `ScriptDemo.m`.
- [x] Gerador logístico do Appendix C (Ives & Helmus 2011), com os sete cenários de
      [E] e [S]: `simulate_logistic()` e `logistic_gradients()`.
- [x] Comparar com o `mcomsimr` nos mesmos parâmetros, levando em conta as diferenças
      de convenção (`sigma`, `alpha[j, k]`). Em `data-raw/validation_mcomsimr.R`: 16
      cenários (dispersão, competição, paisagem), 25 sítios, 12 espécies, 300 passos,
      ambiente constante e 200 réplicas por modelo. Nenhuma das 80 comparações de
      distribuição foi significativa (KS, FDR 5%), e as diferenças nas médias de alfa,
      gama e abundância total ficaram abaixo de 1%. Com extirpação de 0,02, nenhuma das 80
      foi significativa, mas alfa ficou 8% abaixo na grade com dispersão 0,001: nossos
      imigrantes são Poisson, como no código Julia dos autores, e o mcomsimr conserva os
      dispersores (alfa 1,75, contra 1,85 com multinomial e 1,85 no mcomsimr). Falta um
      ambiente que mude no tempo.
- [x] Reproduzir os padrões de Thompson et al. (2020), em
      `data-raw/validation_thompson.R` e `validation_thompson_cc.R`: 100 manchas, 50
      espécies, 2200 passos, 30 paisagens, 15 taxas de dispersão, 6 das 13 larguras de
      nicho, 4 tipos de competição. Os scripts seguem o código dos autores
      (Meta_com_framework v1.0, Zenodo 10.5281/zenodo.3708350), que difere do artigo em
      pontos que mudam o resultado: beta é gama/alfa, e não o gama − alfa da Tabela 1;
      manchas e tempos vazios não entram em alfa; o ambiente é convertido em postos,
      portanto uniforme em [0, 1]; dispersão começa em 0,00001, não 0,0001; a
      competição-colonização reduz os imigrantes à metade. Com isso, em dispersão 1e-5,
      stabilising dá alfa 4,6, gama 42 e beta 9,1 (σ = 0,5) e 6,6, 46 e 7,0 (σ = 10),
      próximos das Figuras 3c e 3d lidas a olho. H1, H2, H3, H5 e H6 aparecem.
- [x] H4 resolvido. A diferença anterior vinha do ambiente, não da definição de beta,
      que já era gama/alfa. A parte da competição aparece (beta espacial em 1e-5, σ = 0,5:
      equal 38,8, mixed 21,5, stabilising 9,1). A parte do nicho depende do cenário, como
      na Figura 4 deles: em stabilising, beta tem pico em σ = 0,215 (14,4 contra 7,0 em
      σ = 10); em mixed o pico é fraco, em σ = 0,5; em equal e dispersão baixa, beta cresce
      com o nicho (16,5, 38,8 e 43,1 em σ = 0,215, 0,5 e 10), por deriva. Nichos estreitos
      só elevam beta em equal a partir de dispersão perto de 0,02. O texto de H4 generaliza
      mais do que a figura mostra. A comparação com as figuras é visual, porque os valores
      publicados não estão em tabela.
- [x] `extirpation` em `simulate_metacommunity()` e `mc_dynamics()`. O cenário
      competição-colonização lembra as Figuras 3g e 3h: com σ = 0,5, gama sobe de 3,6 para
      cerca de 6 perto de dispersão 5e-4; com σ = 10, colapsa em uma espécie perto de 0,02.
- [x] Testes de propriedades com resultado analítico, em `test-analytic.R`: equilíbrios
      com uma e várias espécies, exclusão, fonte e sumidouro, fluxo pelo kernel.
- [x] Artigo "Validating the model of Thompson et al. (2020)" no site, com as divergências
      entre artigo e código, também resumidas no help de `simulate_metacommunity()`.

**Concluído quando:** houver um script em `data-raw/` que reproduza a comparação.

## Métodos de análise

- [x] `pcc_test()` e `phylo_eigenvectors()`: PCC de Peres-Neto & Kembel (2015), vindo do
      pcps-sim, com permutações restritas opcionais. Em `data-raw/validation_pcc.R`, as
      taxas de neither e environment batem com a Tabela 1 (erro tipo I do teste de
      espécies 0,049 e 0,040 contra 0,054 e 0,052). O poder (0,87 e 0,85 contra 0,945 e
      0,947) e o teste de espécies em phylogeny (0,20 e 0,19 contra 0,27 e 0,28) ficam
      abaixo, como esperado sem a árvore deles.
- [x] `pcps_test()`: PCPS de Duarte et al. (2016), com teste de toda a dissimilaridade e de
      cada eixo e `site_permutations`, que o pacote PCPS não tem. `torus_permutations()`
      aceita grade retangular, `symmetries` e `n_perm`. Estatísticas iguais às do PCPS 1.0.9
      (diferença 8e-15 em 20 comunidades) e p compatíveis. Em `data-raw/validation_pcps.R`,
      grade 10 x 10 sem relação real: livre rejeita 29-80% (`p_max` 3-24%), toro 4-7% com
      campos periódicos (`p_max` 1-4%) e 6-25% com bordas (`p_max` 1-10%). O PR ao PCPS
      fica para depois do aceite do manuscrito.
- [x] `phylo_beta_test()` e `phylo_beta()`: framework de beta diversidade filogenética do
      mesmo artigo (db-RDA sobre MPD entre comunidades), que não existia no pcps-sim. F e
      R2 iguais aos do `vegan::adonis2()`, e métrica igual a `picante::comdist()`. Em
      `data-raw/validation_phylo_beta.R`, com 1000 réplicas e o PCC nos mesmos dados, as
      taxas batem com a Tabela 1: neither .005 e .010 (.007, .008), environment .050 e
      .062 (.055, .050), phylogeny .033 e .056 (.035, .039), teste de espécies em
      phylogeny .16 e .19 (.167, .177). Poder .893 e .908 (.937, .939), R2 ajustado .083
      e .060 (.091, .064). O PCC teve poder de .888 e .856 nos mesmos dados, ao contrário
      do artigo, onde fica à frente (.945, .947), e um erro de .058 onde a outra rodada
      deu .068, o que mostra o ruído desses valores.
- [x] Heatmaps do artigo para visualizar os eixos de `scores`: `plot_phylo_beta()`, com `axis`.
- [x] Registrado em `?phylo_beta`: a matriz de MPD não é sempre euclidiana, ao contrário do
      que o artigo afirma, e o pseudo-F não exige isso. A nota antiga de que o picante não traz
      as funções que o artigo anuncia não diz quais são: fica fora do texto até ser conferida.

## v0.3.0: fluxo de estudos de simulação (concluída em 2026-11-07)

Objetivo: reduzir o código que o usuário repete. As vignettes precisaram de laços e
funções escritos à mão, o que indica o que falta no pacote.

- [x] `run_scenarios()`: grades de cenários com réplicas pareadas por árvore e
      paisagem e semente por célula, com paralelismo por `future.apply` em `Suggests`.
- [x] `equilibrium_check()`: diagnóstico de equilíbrio. Com competição filogenética e dispersão alta, a
      exclusão regional continuou até 3000 passos, e hoje nada avisa o usuário.
- [x] Saída em formato longo a partir de `history`: `as.data.frame()`, e `as_tibble()` por cima.
- [x] Benchmark de `mc_dynamics()`, em `data-raw/benchmark_dynamics.R`: 0,7 s para 50
      espécies, 100 sítios e 2200 passos, ou cerca de 19 h por 100 mil simulações em um núcleo.
      O sorteio de números aleatórios (`rmultinom`, `rpois`, `rbinom`) toma cerca de 60% do
      tempo e já é C, então código compilado não traria ganho que justifique. O ABC paraleliza
      por simulação.
- [x] `summary_stats()`: as mesmas estatísticas para dados simulados e observados
      (alfa, beta, gama, ocupação, SAD, equitabilidade, SES MPD e MNTD, decaimento com
      a distância, relação entre traits e ambiente) e, com mais de um tempo, variação
      temporal, como o CV das abundâncias. Em Guzman et al. (2022), estatísticas
      temporais aumentaram o poder explicativo em até 59%, e a dispersão dependeu delas.
- [x] Piloto de identificabilidade, em `data-raw/pilot_identifiability.R`: 30 espécies, grade
      5 x 5, árvore e paisagem fixas, sem erro de detecção; critério combinado antes de rodar,
      R2 >= 0,5 e cobertura de 85% a 95%. Com 5000 simulações, nenhum parâmetro passou: R2 de
      0,8 a 0,98 em dispersão, competição e nicho, mas cobertura de 40% a 90%, que sobe com a
      tolerância. Na sensibilidade fora do plano, com 50 mil simulações e 500 aceitas (tolerância
      0,01), os três passam com um levantamento e com série: R2 0,98 e 0,99 (dispersão), 0,88 e
      0,90 (competição), 0,91 e 0,93 (nicho), cobertura 88% a 92%. A série só estreita o
      intervalo, de 0,13 para 0,10 da amplitude na dispersão. `lambda` não é recuperado
      (R2 <= 0,07). Decisão para o ABC: seguir, com pelo menos 500 aceitas, e dar a `lambda`
      uma estatística de sinal dos traits, que `summary_stats()` ainda não tem. Falta testar com
      árvore e paisagem variando, detecção imperfeita e modelo errado.

**Concluído quando:** as vignettes puderem ser reescritas com essas funções, com bem
menos código. A `phylogenetic-structure` mostra o desenho com `run_scenarios()`; os
resultados salvos não foram refeitos. As vignettes de avaliação de nulos usam laços de outro
tipo e não mudaram.

## v0.4.0: base do hub (concluída em 2026-11-07)

Objetivo: uma interface comum, para que simulação e métodos troquem dados entre si e
com outros pacotes. As versões seguintes dependem desta.

- [x] Estrutura de dados comum: `mc_data()` e `as_mc_data()`, com comunidade, árvore, traits,
      ambiente, coordenadas, áreas e série temporal, e checagem de nomes. As análises aceitam o
      objeto no primeiro argumento, e as simulações guardam a árvore.
- [x] Árvore e distribuições como entrada: `pool` (onde cada espécie pode ocorrer) e
      `initial` (estado inicial) em `simulate_metacommunity()` e `mc_dynamics()`; a árvore
      já entrava por `species_traits()`. Cobre gen3sis, modelos de SSE e reconstruções
      biogeográficas, sem leitor próprio: o formato de saída do gen3sis não foi conferido.
      Pool que muda no tempo (expansões de área) fica para quando houver uso.
- [x] `species_traits()` aceita traits vindos de fora (`optima`, `sigma`, `r_max`), por
      exemplo do mvMORPH; os demais traits saem iguais aos sorteados sem eles.
- [x] Registrar a história dos traits ao longo dos ramos: valores nos nós em `trait_bm()`,
      `species_traits()` (`nodes = TRUE`, sem mudar os outros valores com a mesma semente) e
      `traits_ou()`; `trait_mk()` para trait discreto (Mk), no formato simmap do phytools
      (conferido contra `phytools::sim.history()`); `tip_history_metrics()` com as métricas
      de Luza et al. (2021) como o código de `calc_tip_based_trait_evo()` as calcula (nós,
      não mudanças dentro dos ramos), mais `changes` e `time_since_change` em tempo contínuo.
      O Herodotools não foi instalado (depende do daee, fora do CRAN): a equivalência vem da
      leitura do código, não de uma comparação rodada. Taxas variando entre ramos (relaxed
      BM, BAMM) como modelo gerador ficam para depois.
- [x] Integração com o Herodotools (Nakamura et al. 2024): `as_herodotools()` (entradas com os
      nomes dos argumentos dele e as áreas ancestrais verdadeiras), `region_pool()` e
      `arrival_ages()`, artigo "Known biogeographic histories for Herodotools" e
      `data-raw/validation_herodotools.R` (Herodotools 2.0.1, instalado à parte: o rexpokit do
      BioGeoBEARS não compila neste macOS, usado o binário do CRAN). Resultados: com a história
      verdadeira no lugar dos mapas, as métricas por ponta são idênticas às de
      `tip_history_metrics()` em 30 de 30 histórias. Das pontas, a correlação de Spearman com a
      verdade cai com a taxa: estase 0,97 a 0,69, transições 0,68 a 0,38; as transições nos nós
      perdem as mudanças dentro dos ramos (1,81 contra 3,01 por ponta na taxa mais alta). Idade
      de chegada com áreas verdadeiras: 0,70 a 0,89; com áreas por ML (`ape::ace`): 0,24 a 0,58,
      porque um nó profundo errado muda a idade de um clado inteiro. Herodotools fora do
      `Suggests` por não estar no CRAN.

**Concluído quando:** um cenário simulado passar pelo Herodotools e recuperar a
história conhecida.

## v0.5.0: métodos de metacommunity phylogenetics

Objetivo: levar para o R os métodos que só existem em MATLAB. O código está nos
suplementos de Peres-Neto et al. (2012) e Leibold et al. (2010), e nenhum dos dois
métodos está em pacote de R.

Vem antes do ABC por decisão de 2026-11-06. O código está em
`~/Downloads/ECOL_93_175/` (Peres-Neto et al. 2012) e em
`~/Downloads/ele_1523_sm_datas2.doc` (Leibold et al. 2010), fora do repositório porque
os suplementos não trazem licença. As comparações usaram Octave 11.3, instalado nesta
máquina.

- [x] Peres-Neto et al. (2012): `phylo_vector()` e `trait_partition()`, conferidas no Octave e nas Tabelas 1
      a 4 do artigo.
- [x] Tabela 5 com `logistic_gradients()` (Apêndice C), em `data-raw/validation_table5.R`.
      Frações com os vetores que geraram os dados e rejeição de E e S pelo teste global de
      `phylo_vector()`: o padrão entre cenários se repete, e I, III e IV ficam mais perto. A
      razão de variância não bate (0,62 contra 0,39 a 0,94), e o cenário VII tem `b` maior.
      A diferença está no artigo de validação; para fechar, só perguntando aos autores.
- [x] Semiparciais de `trait_partition()`, do Apêndice A, com testes. O suplemento não
      tem código para elas, então a conferência foi feita contra regressões separadas.
- [x] `node_test()` roda nós em paralelo com `cores`. Medido em um núcleo, com 100 sítios e
      999 permutações: 0,5 min com 100 espécies, 1,7 com 300 e 3,9 com 600, linear nos nós.
      Não é lento, então `node_fit()` não foi vetorizado.
- [x] Leibold et al. (2010): `node_test()` e `node_contrasts()`. R2, R2 ajustado e F iguais
      ao MATLAB no Octave 11.3 (diferença abaixo de 1e-15), e os 13 contrastes de
      `CreateNodeContrasts` iguais. Em `data-raw/validation_node.R` (refeito na v0.5.2, com F), por nó, sites .053 e
      espécies .049 sob a nula dupla; p máximo .049 só com ambiente e .021 só com filogenia.
      Poder por nó .07, e ao menos um nó em 76% das réplicas, mas 24% sob a nula: é preciso
      corrigir para testes múltiplos. Os erros do original foram confirmados na leitura e
      corrigidos, com a diferença descrita em `?node_test`.
- [x] Semiparcial, FDR e eventos biogeográficos conforme o texto de Leibold et al. (2010).
      O MATLAB multiplica pelo módulo da inclinação e usa preditores brutos; o artigo
      divide e usa padronizados, e seguimos o artigo.
- [x] Codificação com comprimentos de ramo, na seção III do suplemento de Leibold et al.:
      `branch_lengths = TRUE`, conferida na Figura S2.
- [x] `plot_node_test()`: Figura 3 de Leibold et al. (2010), em filograma ou, para árvores
      muito grandes, leque, com os nomes das variáveis ao longo dos raios, com
      `ape::plot.phylo()`, sem dependência nova. O artigo não define a família do FDR:
      corrigimos preditores e eventos separadamente, o que está no help.
- [x] `plot_trait_partition()`: frações em barras, com negativas como zero (Legendre 2008).
      A tabela mantém os valores negativos.
- [x] `node_test()`: refúgios, 9999 permutações por padrão, permutação conjunta do artigo
      (`p_joint`, que rejeita 76% dos nós só com ambiente, contra 4,8% de `p_max`).
- [x] `daphniids`: dados de Leibold et al. (2010), refeitos da pesquisa pública da EPA.
- [x] Conferência de espécies por nome nas funções com árvore, traits ou contrastes.
- [x] Saíram `fourth_corner_test()` (fica no `ade4`) e `spatial_eigenvectors()` (fica no
      `adespatial`). `trait_partition()` cita ter Braak (2019) e o TraitEnvMLMWA para o
      teste de cada traço.
- [x] `plot_node_heatmap()`: nós por preditores com `semipartial`, como a Tabela S1.
- [x] Estilo tidyverse na entrada e na saída, com álgebra de matrizes por dentro: tibbles
      numéricos entram, `as.data.frame()` e `as_tibble()` saem. Permutações em paralelo com
      `future.apply`, em `Suggests`, sob o `plan()` do usuário; `node_test()` segue com `cores`.
- [x] Testar a tradução contra o MATLAB em dados fixos, antes de corrigir qualquer
      coisa: R e Octave concordam até 1e-13, e as divergências estão nos helps.
- [x] Artigo "Analysing clades and traits in a metacommunity", com dafinídeos e simulação.

Situação em 2026-11-06: marco fechado na 0.5.1. Sobra o estilo tidyverse e o Futureverse,
se ainda fizerem sentido.

**Concluído quando:** a tradução reproduzir os resultados do MATLAB e os testes
reproduzirem o erro tipo I publicado.

## v0.5.2: biogeografia em `node_test()` e auditoria

Objetivo: tirar a maior barreira de uso do método de Leibold et al. (2010), que é montar a
matriz de eventos, e revisar o pacote do ponto de vista de quem usa. Vem antes de v0.3.0 e
v0.6.0 por decisão de 2026-11-06.

- [x] `biogeography_matrix()`: matriz de eventos a partir da área de cada site ou de uma
      matriz 0/1 de sites por áreas. Pares de áreas definidos a priori são o caminho
      principal, como em Leibold et al. (2010) e Leibold & Chase (2018, cap. 9), que usam
      hipóteses publicadas. "Cada área contra as outras" fica como opção exploratória, com
      fora das áreas como 0 (script de 2022) ou como -1 (script de 2015), e aviso quando um
      lado tiver poucos sites. Os dois scripts estão em `legacy/` do af-frogs-climate-rarity.
      Sobrepor polígonos e grade fica para o usuário, com `sf::st_join()`, sem dependência
      nova. Reproduz os seis eventos de `daphniids` e a tabela de refúgios de 2022 (com 5 como 1).
- [x] `site_permutations` em `node_test()`: para eventos, cada linha move as áreas pelos
      sites. Em `data-raw/validation_events.R` (1000 réplicas, toro de 15 por 15, blocos de 25
      sites sem efeito, comunidade seguindo ambiente oculto autocorrelacionado), rejeição por
      nó a 5%: .451 com permutação livre e .048 com o toro; ao menos um nó em 89% e 26% das
      réplicas, sem correção. Só testado com toro.
- [x] Nulo restrito para paisagens com borda: `patch_permutations()` e `area_permutations`
      em `node_test()`. Nenhum método publicado cobre permutação exata, bordas, área
      irregular e variável binária, então é método próprio. Em
      `data-raw/validation_edges.R` (200 réplicas; quadrado, faixa e banda irregular; blocos,
      áreas em manchas e área contra as outras), rejeição por nó a 5% sem efeito: .030 a .063,
      contra .44 a .55 com permutação livre. Descartados no protótipo: curva de Hilbert com
      deslocamento cíclico (.16 a .19), rotação rígida (lenta; .06 a .12) e mancha única por
      lado (conservador em áreas com várias manchas, .02 a .03, e menos poder). Falta aplicar a
      dados reais.
- [x] `plot_node_test()`: triângulos para eventos, como a Fig. 9.5 de Leibold & Chase
      (2018), além da cor.
- [x] Citar Leibold & Chase (2018) em `biogeography_matrix()`.
- [x] Auditoria do código com foco em experiência de uso, com dois revisores (2026-11-06).
      Corrigidos os erros que mudam resultados e as falhas silenciosas: F nos sorteios de
      `node_test()` (a estatística parcial do MATLAB dependia da escala), nós com clado
      ausente, espécies ausentes em `pcc_test()`, `phylo_beta_test()` e `trait_partition()`,
      `cores`, sites por nome, checagens em `R/checks.R` e `permutations` nos resultados.
- [x] `plot_trait_partition()`: barras sobrepostas, como a Fig. 1 de Desdevises et al. (2003).
- [x] Exemplos, README e vinhetas com `treeshapesim::sim_tree_shape()`. Testes seguem com
      `ape`. Os `.rds` das vinhetas foram refeitos: com a árvore nova, SES MPD fica em -1,3
      (era -1,7), e o toro rejeitou .08 e .07 (eram .035 e .050, com erro padrão .015).
- [x] Nota "Adapted from the MATLAB" em `node_test()`, `node_contrasts()`, `phylo_vector()` e
      `trait_partition()`. Ver a pendência sobre licença do código.
- [x] Barra de progresso com `purrr` e estimativa de tempo em sessões interativas, em
      `node_test()`, `pcc_test()`, `phylo_beta_test()`, `phylo_vector()` e `trait_partition()`.
      `purrr` entra em `Imports`.

Da auditoria, para antes da 1.0:

- [x] Descoberta: tabela de pergunta para função no README, página `?comdynamite`, novo Title
      e Description, artigo sobre preparar dados próprios e artigo de `pcc_test()` e
      `phylo_beta_test()`.
- [x] Nomes: `level` no lugar de `alpha` para significância, `n_perm` com padrão 999, `env` em
      `phylo_vector()`, `p_species`, `p_partial_max`, `r2_mean` em `pcc_test()` e `r2_adjusted`
      em `phylo_beta_test()`. `statistic` fica: é a estatística de cada teste, descrita na ajuda.
- [x] Métodos `summary()`, `plot()` e `as.data.frame()` para os resultados, em `R/methods.R`.
- [x] R2 de `node_test()` acima de 1 (herdado do original): aviso na ajuda e `r2_weighted`.
- [x] `site_permutations` em `trait_partition()`, e sites vazios com permutações restritas em
      `pcc_test()`, `phylo_beta_test()`, `phylo_vector()` e `trait_partition()`.
- [x] Baixa gravidade: ramos zero em `node_contrasts()`, nomes repetidos e heatmap só com NA,
      `logistic_gradients()` com poucos sites, `mc_dynamics()` com um site, `r2` NaN em
      `pcc_test()`, `n_draws`, erro dentro de `cores`, checagens do simulador.
- [x] Formas dos traits, na auditoria com Fable (2026-11-07): `alpha` alinhado por nome em
      `simulate_metacommunity()`, `as_traits()` para matriz e tibble, aviso de nomes inventados.

Situação em 2026-11-06: marco fechado na 0.5.2. O nulo com borda e a descoberta saíram na
0.6.0.

Opção, sem data: calanoides como segundo exemplo, se as espécies estiverem nas tabelas
públicas da EMAP. A topologia sairia da Fig. 9.5b, e os resultados do livro não foram
publicados em artigo, então seria exemplo novo, sem afirmar que reproduz o livro. O texto
do livro diz que os dafinídeos não mostram restrição biogeográfica, mas a Fig. 9.5a tem
vários triângulos: vale perguntar aos autores.

**Concluído quando:** a matriz dos refúgios e a dos dafinídeos saírem da função, o erro tipo
I dos eventos estiver medido com sites agregados no espaço e os achados da auditoria
estiverem resolvidos ou registrados aqui.

## v0.6.0: nulo com borda e descoberta (concluída em 2026-11-06)

`patch_permutations()`, `area_permutations` em `node_test()`, novo título e documentação por
pergunta. Os detalhes estão na seção da v0.5.2, onde os itens nasceram.

- [x] `patch_permutations()` com autocorrelação na escala da paisagem (bandas de 6, 12 e 25
      sites, 200 réplicas, sem efeito): rejeição de .022 a .060, contra .49 a .75 com
      permutação livre e .056 a .095 com o toro no quadrado. Conservador para área contra as
      outras (.022 a .039). Em `data-raw/validation_edges_scale.rds`.
- [x] Aplicado aos refúgios da Mata Atlântica, fora deste repositório (dados não publicados).
      Com células esparsas, o resultado depende muito de `gap`.
- [x] `patch_permutations()` não preserva a distância entre os lados do evento: resolvido em
      parte com `rigid_permutations()` (movimento rígido do evento inteiro; `independent` para
      lados distantes). Lados próximos: .022 a .048. Eventos que ocupam a extensão da paisagem
      ficam sem layouts ou, com `independent`, liberais (.06 a .09; .12 com autocorrelação longa).
      Não há nulo válido para lados distantes sem conhecer o alcance da autocorrelação.
      Em `data-raw/validation_rigid.R`. Registro do problema: com áreas fixas
      no quadrado 20 por 20: adjacentes .004 a .02; cantos opostos e metades .06 a .10 (.065 e
      .064 na minha conferência, 150 réplicas; .078 a .101 na auditoria, 60 réplicas). A
      validação sorteava as áreas, o que tira a média disso. Candidato: mover o conjunto de
      manchas de um evento junto, com deslocamento e rotação comuns.
- [x] Ramo fora de grade de `rigid_permutations()`: agora enumera rotações (passo que move o
      sítio mais distante em cerca de `gap`, de 4 a 360), reflexões e translações para cada
      sítio, sobre um raster com fechamento (buracos contam como dentro, a borda não), em vez
      de sortear movimentos; o conjunto de encaixes é fixo. Em `data-raw/validation_offgrid.R`
      (150 réplicas, 99 permutações): áreas adjacentes .047 e .059; dois blocos aleatórios em
      quatro paisagens .068 a .080, liberal; `patch_permutations()` .012 a .063; poder 23%
      contra 12%. Fica pendente: blocos aleatórios fora de grade acima de .07, e eventos
      grandes em grade com jitter são lentos (225 sítios em 750: 701 encaixes em 12 min).
      Nos refúgios, a parte grande tem um único encaixe, de forma determinística. Conferência
      independente (região em L, `data-raw/validation_nulls_check.R`): .070 e .087 fora de
      grade, contra .055 e .066 com `patch_permutations()`; grade esparsa .057 e .061.
- [x] Critério para `gap`: mediana da distância ao quarto sítio mais próximo (passo da grade
      cheia; espaçamento com sítios esparsos). `patch_permutations()` usa 1,5 vez isso (cantos
      com 25% a 40% da grade: .056 a .059, contra .15 a .17 com o padrão antigo; holdout .021
      a .059). `rigid_permutations()` usa o valor direto e não erra mais com poucos buracos
      (.055 a .061 com 40% a 80% dos sítios). Em `data-raw/validation_gap.R`. Cantos opostos
      esparsos e clusters seguem sem encaixes para o rígido.
- [x] Testes das frações `a` e `c` de `trait_partition()`. O MATLAB permuta ambiente e espaço
      com a mesma ordem de sítios e compara a própria fração (calcula o F parcial e não usa).
      Em vez de trocar a permutação, trocamos a estatística pela parcela do componente que o
      outro conjunto deixa (monótona no F parcial), o que mantém `site_permutations` e a
      permutação de espécies. No tuning (200 réplicas, grade 8 por 8) ela equivale ao
      Freedman-Lane e à permutação só de `env`; resultados do holdout (500 réplicas) em
      `data-raw/validation_fractions.R` e no artigo de validação. Limites registrados: com
      `env` estruturado em escalas fora de `space` o teste livre de `a` infla, e o toro
      segura; o toro é quase degenerado para `c`, porque mover os MEMs rigidamente quase não
      muda o espectro deles.
- [x] Teste da fração compartilhada `b` com `shared_null`: `"msr"` (Moran spectral
      randomization de `env` a partir de `coords`, via `adespatial::msr()` com a base de
      `dbmem()`) ou `"sites"` (`env` sozinho movido por `site_permutations`). Só testa se `env`
      estiver relacionado a `space` (`level`), como `envspace.test()`. p de sítios do MSR, p de
      espécies da permutação de traits já existente, `p_max` o maior; a permutação de espécies
      responde se o trait importa, não se a estrutura compartilhada é espúria, e o teste de
      `p_max` precisa dos dois (com trait oculto, o MSR sozinho rejeita 54%). Holdout com
      `env` estruturado e sem efeito: permutação livre de `env` .48, MSR .056, toro .11
      (liberal, sem explicação; fica como limite). Poder: 1.0 só ambiente, .43 com espaço
      agindo também, .28 com ambiente fraco. Referência verificada: Bauman et al. (2019,
      Oikos); Clappe et al. (2018) trata da partição com MSR, não do teste de `b`.
- [x] Nulo restrito para `c` em `trait_partition()`: `site_permutations` move `env` e `space`
      juntos, e mover MEMs de uma grade quase não muda o que explicam; no holdout o toro foi
      liberal para `c` com `env` forte. Candidato: mover só `env` nos testes de `a` e manter
      `c` com permutação livre, ou documentar `c` como descritivo com nulo restrito.
      Resultado (`data-raw/validation_fractions_null.R`): cada conjunto passa a se mover
      sozinho (`env` para `a`, `space` para `c`), o que segura o nível nos dois geradores
      (`a` .037 a .071 e `c` .030 a .080, contra .20 e .24 da permutação livre) e derruba a
      premissa de degenerescência (os MEMs movidos têm em média 30% da variância fora do span
      dos originais), mas o teste restrito de `c` só tem poder com campo no span de `space`
      (1.0; .25 a .28 com campo exponencial que `space` segura em parte, contra .98 livre),
      porque pergunta se esses mapas explicam mais do que os mesmos mapas em outro lugar da
      grade, não se o resíduo tem estrutura espacial.

## v0.6.3: base do hub (concluída em 2026-11-07)

Publica o marco v0.4.0: `mc_data()` já saíra na 0.6.2; entram traits vindos de fora em
`species_traits()`, `pool` e `initial` nas simulações, valores nos nós, `trait_mk()`,
`tip_history_metrics()` e a integração com o Herodotools (`as_herodotools()`,
`region_pool()`, `arrival_ages()`), além de `trait_k` em `summary_stats()` (marco do ABC). DOI da versão:
`10.5281/zenodo.24046171`.

## v0.6.2: fluxo de simulação, PCPS e estrutura comum (concluída em 2026-11-07)

Publica o marco v0.3.0 (`run_scenarios()`, `summary_stats()`, `equilibrium_check()`, saída
longa, benchmark e piloto de identificabilidade), `pcps_test()` com o toro retangular, as
pendências das etapas anteriores, o nulo restrito por conjunto em `trait_partition()` e o
primeiro item da v0.4.0, `mc_data()`. DOI da versão: `10.5281/zenodo.24044811`.

## v0.6.1: nulo rígido e testes das frações (concluída em 2026-11-07)

`rigid_permutations()`, enumeração fora de grade, critério para `gap`, testes das frações
`a`, `b` e `c` de `trait_partition()`, `plot_phylo_beta()`, métodos de resultados e os nomes
padronizados. Substitui a `v0.6.0` para citação: a tag `v0.6.0` traz as taxas antigas de
`patch_permutations()` e não teve release. O nulo restrito para `c` fechou depois, na
versão em desenvolvimento (item da v0.6.0).

## Inferência por ABC (planejada como v0.6.0, número a definir)

Objetivo: estimar parâmetros e comparar modelos a partir de dados observados. Depende do
piloto da v0.3.0, já feito, e da estrutura comum da v0.4.0.

- [x] Estatística de sinal filogenético dos traits em `summary_stats()`: `trait_k`, K de
      Blomberg, igual ao `picante::Kcalc()`. Piloto repetido (`pilot_identifiability_k*.rds`):
      com 5000 simulações, nada passa de novo (cobertura baixa); com 50 mil e tolerância 0,01,
      `lambda` passa por pouco, R2 0,53 e cobertura 90-91%, com intervalo de 61-63% da amplitude.
      Os outros três seguem como antes. Ressalva: o K é dos ótimos, que não se observam no
      campo; com um trait medido, a recuperação de `lambda` deve ser pior. Com 30 espécies o K
      é ruidoso; árvores maiores devem ajudar (não testado).

Ler antes de implementar: Guzman et al. (2022), lido em parte; mcfly (Duarte et al. 2023); CAMI (Ruffley et al.
2019), pacote de R com random forests e ABC para montagem de comunidades; van der Plas
et al. (2015). Definir o que o comdynamite acrescenta a esses trabalhos.

Pontarp et al. (2019), lido com o Appendix S4: propõem um ciclo
de construir o modelo, analisá-lo sem dados, estimar parâmetros e avaliar. O caso de
estudo usa ABC por rejeição com 14 mil simulações por modelo, distância euclidiana entre
estatísticas (riqueza, abundância, distribuição de traits, NRI, NTI) escaladas pelo desvio
padrão, priors amplos e dois modelos, escolhidos pela razão de aceitação, tudo em dados
simulados e só com a comunidade final. O limiar sai do gráfico das distâncias e de uma
análise de sensibilidade; vale copiar esse passo, mas ele muda entre os trechos (1 na
Fig. 5, 3 e 5 no apêndice). O que levar para cá:

- Se um processo não é recuperado em dados simulados, não deve ser inferido em dados
  reais. É o papel do piloto da v0.3.0 e da recuperação abaixo.
- Traits e filogenia juntos separaram processos com sinais parecidos, então
  `summary_stats()` deve cobrir os dois.
- O modelo errado deu estimativas aparentemente boas, com a largura de nicho
  superestimada. Uma boa estimativa não valida o modelo: daí os PPCs.
- Eles não selecionaram estatísticas e citam Prangle et al. (2014) para isso. A razão de
  aceitação depende do limiar e das estatísticas (Robert et al. 2011), então preferir
  `abcrf` e medir o erro de classificação em dados simulados.
- O que o comdynamite acrescenta: espaço com muitos sítios, pacote em R (o código deles
  é MATLAB), cobertura dos posteriores e checagem do modelo, que o caso de estudo não faz.

- [ ] Tabela de referência: sortear parâmetros dos priors e simular, com a árvore, o
      ambiente e as coordenadas observados e com detecção imperfeita.
- [ ] Estimação e escolha de modelos com `abc` (Csilléry et al. 2012) e `abcrf` (Pudlo
      et al. 2015), em `Suggests`. Random forests lidam bem com muitas estatísticas
      correlacionadas e ordenam a importância delas.
- [x] Custo, medido em `data-raw/benchmark_dynamics.R`: 0,7 s para 50 espécies, 100 sítios
      e 2200 passos, como em Guzman et al. (2022), ou cerca de 19 h por 100 mil simulações em
      um núcleo; a grade 5 x 5 do piloto levou 1,1 min por 5200 simulações em 6 núcleos.
- [ ] Recuperação em dados simulados com parâmetros conhecidos: viés e cobertura.

Avaliação do modelo, sempre em conjunto com a inferência:

- [ ] Posterior predictive checks, com estatísticas que não entraram no ajuste.
- [ ] LOO sobre a tabela de referência: `cv4abc()` e `cv4postpr()` do `abc`, com cada
      simulação como pseudo-observação; `gfit()` para adequação global.
- [ ] K-fold por sítios ou blocos espaciais: ajustar com parte dos sítios e prever as
      estatísticas dos demais. Blocos evitam que a autocorrelação espacial infle o acerto.
- [ ] Diagnóstico de Pareto k dos pesos de importância, com `loo::psis()`, para saber
      quando a posterior ponderada depende de poucas simulações.
- [ ] PSIS-LOO padrão (Vehtari et al. 2017) precisa de log-verossimilhança por observação,
      que o ABC não dá. Avaliar likelihood sintética (Wood 2010) por sítio ou por
      estatística, que daria esse termo de forma aproximada.

Limites a deixar claros: um bom ajuste não prova o mecanismo, e dados de um único
levantamento podem não informar a dispersão.

**Concluído quando:** uma vignette levar dados reais do ajuste à avaliação, com
cobertura medida em dados simulados.

## v0.7.0: evolução de traits e modelos nulos

Objetivo: ampliar os modelos de traits e os nulos com dinâmica histórica, chamando os
pacotes de referência em vez de reimplementá-los.

- [ ] Mais opções de evolução de traits, multivariadas e com covariância entre traits,
      pelo mvMORPH (Clavel et al. 2015).
- [x] Taxas que variam entre ramos como modelo gerador. Avaliação da literatura (2026-11-07):
      um gerador com a verdade conhecida, não um método de inferência; RRphylo, phytools,
      motmot, bifrost e corHMM são GPL e só entram em `Suggests` ou `data-raw`. O RRphylo não
      tem modelo discreto (aceita estados codificados como números, tratados como contínuos).
      BayesTraits e RevBayes adiados (binários externos; licença do BayesTraits não
      declarada); BAMM sem manutenção desde 2017.
      `trait_rates()`: modos `"shift"` (mudanças herdadas ou só no ramo), `"gbm"` (log da taxa
      por BM) e `"iid"` (lognormal por ramo), com os valores nos nós e as taxas por ramo.
      `trait_mk()`: `alpha`/`edge_rates` (Γ por ramo, como `fitgammaMk`) e `rate_classes`
      (classes ocultas, como corHMM). Validação em `data-raw/validation_rates.R`.
      `shift_rate` (mudanças por processo de Poisson ao longo dos ramos, como no prior do
      BAMM) e `relative = FALSE` (cada regime com taxa própria, como no BAMM e no bifrost).
      Validação com o BAMM 2.5.0 em `data-raw/validation_bamm.R` (binário local, GPL, só
      executado): Spearman 0,51 entre taxa média posterior e verdadeira por ramo; probabilidade
      marginal de mudança 0,18 nos ramos com mudança verdadeira contra 0,02 nos demais.
- [ ] Artigo (vignette) sobre taxas por ramo. Citar obrigatoriamente Hubert et al. (2015),
      modelos de especiação em metacomunidades, e Januario et al. (2025), persistência
      metapopulacional contra taxas de especiação nas pontas (sinal filogenético, mas sem
      relação com a especiação). Ponte: os dois tratam de taxas de diversificação; aqui as
      taxas são de evolução de traits, e a pergunta análoga é se a dinâmica da metacomunidade
      deixa rastro nas taxas estimadas por ramo. Para o papel da montagem de metacomunidades
      na evolução de traits, citar também Leibold et al. (2022), evolução e montagem de
      comunidades através de escalas espaciais, e McPeek (2017), *Evolutionary Community
      Ecology*. Notas de leitura e cenários de simulação em `dev/notas-artigo-taxas.md`.
      Incluir os resultados da validação e o limite de identificabilidade (uma realização por
      ramo).
- [ ] Ponte Leibold–Januario como eixo do artigo (detalhes em `dev/notas-artigo-taxas.md`).
      Leibold et al. (2022) apontam que os modelos macroecológicos não incorporam quanto a
      adaptação local e a montagem alteram especiação e extinção e, por elas, as filogenias.
      Januario et al. (2025) testaram uma versão disso, a persistência metapopulacional
      contra a taxa de especiação nas pontas, e acharam desacoplamento. Leibold et al. (2019)
      e Urban et al. (2020) dão mecanismos para esse desacoplamento: efeitos de prioridade
      mediados pela evolução tornam as espécies equivalentes na escala regional (neutralidade
      regional) e a evolução em geral amortece os padrões espaciais (85% dos estudos
      revisados por Urban et al.). Hipótese de trabalho (interpretação nossa): a dinâmica
      local deixa pouco rastro nas taxas entre espécies porque a divergência adaptativa fica
      entre populações, não entre espécies. O comdynamite pode testar isso com verdade
      conhecida: metacomunidade → história por linhagem (`history_rates()`) → taxas por ramo e
      traits (`trait_rates()`) → montagem → estimadores de taxa por ramo, com traits de
      indivíduos e populações (estatísticas T) para separar a escala intraespecífica da
      interespecífica. Cenários: frequência de extinção de manchas e dispersão (os dois eixos
      da Fig. 5 de Leibold et al. 2019), autocorrelação espacial do ambiente, ordem de
      chegada, as três taxas relativas da Fig. 1 de Leibold et al. (2022). Cuidado: Leibold et
      al. e Urban et al. falam de escala ecológica; Januario et al. de taxas de especiação, e
      nós de taxas de evolução de traits. O artigo deve dizer isso.
- [x] `food_web_dynamics()`: modelo eco-evolutivo de três níveis tróficos do capítulo 3 de
      McPeek (2017), eq. 3.7–3.17 e 3.20–3.21, reimplementado das equações publicadas (o
      código MATLAB do autor não tem licença declarada e não foi consultado). Abundâncias e
      médias dos traits integradas juntas por RK4 de passo fixo em R base; gradientes de
      seleção por diferenças centrais, o que permite as três formas de coeficiente de ataque
      (gaussiana, logística, Michaelis-Menten) e a omnivoria sem derivadas à mão. Níveis
      ausentes e introduções em tempos dados (`introductions`), para testar a previsão de
      que a evolução dispara após uma invasão e desacelera até o equilíbrio e que o predador
      enfraquece a seleção sobre o recurso; a saída guarda os gradientes e |dz/dt| por nível,
      para virar taxa por linhagem. Validação em `data-raw/validation_food_web.R` (fig. 3.2,
      3.3 e 3.11). Duas leituras da eq. 3.20 documentadas na ajuda: o termo de perda do
      recurso para o predador sem `w`, e `m` como função dos traits do consumidor e do
      predador (box 3.1).
      Falta: traduzir as taxas por nível em `edge_rates` para `trait_rates()`, várias
      espécies por nível (cap. 5) e o cenário 4 das notas.
- [x] Função que traduz a história ao longo da árvore em taxas por ramo: `history_rates()`,
      com taxa por estado (`"state"`, como O'Meara et al. 2006), explosão após cada mudança
      que decai exponencialmente (`"burst"`, cenários 1 e 2 das notas) e dependência do número
      de linhagens no mesmo estado (`"sympatry"`, o DDexp com geografia de Drury et al.
      2016); `trait_rates(edge_rates = )` recebe as taxas. `"sympatry"` validado contra o
      RPANDA em `data-raw/validation_sympatry.R`: matrizes de covariância iguais às do DDexp
      com geografia (diferença máxima 3e-15) e r recuperado por `fit_t_comp()` (média -0,30
      para -0,3 verdadeiro, 0,005 para 0). `CreateGeoObject()` falha em mapas sem mudança de
      região (bug do RPANDA, objeto `old.Dist`); essas histórias são sorteadas de novo.
      Falta ligar a saída de `simulate_metacommunity()` (entradas e saídas de parceiros,
      cenários 5 e 8).
- [ ] Indivíduos a partir das simulações (variação intraespecífica e adaptação local por sítio)
      e estatísticas T de Violle et al. (2012) via `cati::Tstats()` em `Suggests` (cenário 10
      das notas).
- [ ] Funções que usam traits, como `assemble_community()`, com mais de um trait e com
      restrições na estrutura entre traits, como em `mvSIM()` e nos ajustes do mvMORPH:
      sem covariância (`"diagonal"`), taxas iguais (`"equal"`), variâncias ou correlações
      iguais (`"variance"`, `"correlation"`), matrizes proporcionais (`"proportional"`) ou
      compartilhadas (`"shared"`), e matriz definida pelo usuário. Definir como cada processo
      usa vários traits: distância multivariada na similaridade limitante, nicho
      multidimensional no filtro.
- [ ] Modelo nulo dinâmico do DAMOCLES (Pigot & Etienne 2015): colonização e extinção
      local ao longo da árvore.

- [ ] Rodada de UX antes de fechar a 0.7.0 (pedido em 2026-11-08): revisar nomes e
      argumentos das funções novas (`trait_rates()`, `trait_mk()`, `history_rates()`,
      `food_web_dynamics()`) quanto à consistência com o resto do pacote; pensar em funções e
      métodos de plot (taxas por ramo pintadas na árvore, trajetórias de `food_web_dynamics()`,
      ordenação PCPS); melhorar a documentação (exemplos, `@seealso`, uma vignette de taxas
      por ramo) e as mensagens de erro.

**Concluído quando:** cada opção tiver teste e exemplo em vignette.

## Especiação acoplada à metacomunidade (número a definir, depois da v0.7.0)

Objetivo: fechar o limite da ponte Leibold–Januario. Hoje a árvore é dada e a especiação não
depende da metacomunidade. Plano completo, com o modelo, os modos de especiação, a validação e
as fases, em `dev/plano-especiacao.md`.

- [ ] Fase 1: traits por população (espécie × sítio) com seleção local e fluxo gênico em
      `mc_dynamics()`, sem especiação; reproduzir qualitativamente a neutralidade regional de
      Leibold et al. (2019).
- [ ] Fase 2: divergência entre populações e especiação (`"allopatric"`, `"ecological"`,
      `"protracted"`, `"point"`), com registro filogenético (árvore completa e reconstruída) e
      extinção global; validação no limite neutro e contra o PBD.
- [ ] Fase 3: saídas (simmap da distribuição de cada linhagem, métricas de persistência de
      Januario et al. 2025, traits por população) integradas a `history_rates()`,
      `trait_rates()` e às estatísticas T.
- [ ] Fase 4: comparação com o gen3sis (GPL, só em `data-raw`) e medida de custo.
- [ ] Fase 5: experimentos (persistência contra especiação nas pontas com verdade conhecida;
      neutralidade regional contra taxas por ramo) e artigo.
- [ ] Decidir: função nova (`simulate_diversification()`, provisório) ou argumento em
      `simulate_metacommunity()`; modos da primeira versão; número da versão.

**Concluído quando:** as fases 1 a 4 tiverem testes, validação em `data-raw` e uma vignette.

## v0.8.0: traits, filogenia e biogeografia

Objetivo: transformar o desenho de Benício et al. (2021) em um fluxo do pacote e testar,
com a verdade conhecida, se ele separa tempo para especiação, taxas de diversificação
e conservação de nicho. Scripts e dados do artigo estão no Dryad
(https://doi.org/10.5061/dryad.sf7m0cg65).

Simular o cenário:

- [ ] Árvore e regiões geradas juntas, por exemplo com um simulador de GeoSSE, em cenários
      que variem especiação, extinção, dispersão e tempo de colonização por região.
- [ ] Filtragem por traits discretos, como os modos reprodutivos. `trait_mk()` (v0.4.0)
      já simula o trait com a história; falta o modelo de comunidade usá-lo.
- [ ] Duas regiões com pools próprios e dispersão limitada entre elas, usando a
      entrada da v0.4.0.
- [ ] Amostragem imperfeita: sítios, dias de coleta e espécies não detectadas.

Analisar, chamando os pacotes de referência sobre a estrutura comum:

- [ ] Riqueza por rarefação e extrapolação (`iNEXT`).
- [ ] Sinal filogenético em traits discretos: Rao decomposto nos nós (`adiv`, `ade4`).
- [ ] Diversificação por região: GeoHiSSE (Caetano et al. 2018), no lugar do GeoSSE
      usado no artigo, e testes baseados nas pontas, como FiSSE (Rabosky & Goldberg
      2017) e ES-sim (Harvey & Rabosky 2017).
- [ ] Curvas de linhagens no tempo por região (`ape`).
- [ ] DEC (BioGeoBEARS; Matzke 2014) e tempo de colonização por idades crown e stem.
- [ ] MPD, MNTD e análogos funcionais, com o pool definido em três extensões: só uma
      região, só a outra e as duas juntas (`picante`).
- [ ] Relação da diversidade com clima, clima passado e topografia nas três extensões.
- [ ] Diversificação em função de clima e topografia, com o RPANDA (Morlon et al. 2016).

Depois:

- [ ] Avaliar o fluxo: em que cenários ele recupera a causa simulada e com que frequência
      erra. Incluir árvore incompleta e poucos sítios, limitações citadas no artigo, e
      os falsos positivos de modelos SSE (Rabosky & Goldberg 2015). Comparar GeoSSE,
      GeoHiSSE e testes baseados nas pontas nos mesmos cenários.
- [ ] ClaSSE (Goldberg & Igić 2012) e outros modelos de Markov com k estados.
- [ ] SBEARS, do Herodotools. Conferir o que ele faz: entrou aqui como modelo contínuo,
      e o nome é Site Based Estimation of Ancestral Range of Species.

**Concluído quando:** um estudo de simulação mostrar quando o fluxo acerta e quando erra.
Esse estudo pode virar artigo próprio.

## v0.9.0: mais processos

Objetivo: ampliar os cenários do simulador. Escolher pelo uso real.

- [ ] Vários eixos ambientais (hoje há um só).
- [ ] Paisagens irregulares ou lidas de coordenadas, além da grade.
- [ ] Null models restritos com bordas para preditores contínuos e para o simulador. Em
      `pcps_test()`, o toro numa grade com bordas rejeitou 6% a 25% sem relação real.
      `patch_permutations()` cobre só eventos; candidato: Constant Moran permutation
      (Arthur 2024), a testar.
- [ ] Gerador de ambiente temporal fiel ao de Thompson et al. (2020), como alternativa
      a `env_over_time()`.

**Concluído quando:** os itens escolhidos tiverem testes e exemplo em vignette.

## v1.0.0: CRAN e artigo

Objetivo: API estável e pacote no CRAN.

- [ ] Congelar a API e registrar as funções depreciadas no `NEWS.md`.
- [ ] Tirar `Remotes` do `DESCRIPTION`: o CRAN não aceita. `treeshapesim` precisa estar
      no CRAN, ou ficar só em `Additional_repositories` com o CI ajustado.
- [ ] `cran-comments.md`; rodar `R CMD check --as-cran`, `rhub` e `win-builder`.
- [ ] Conferir o tempo dos exemplos e das vignettes.
- [ ] Submeter ao CRAN e criar a tag `v1.0.0`.
- [ ] Artigo de software. Periódico e momento a definir; depende do manuscrito do
      pcps-sim (ver abaixo).

**Concluído quando:** o pacote estiver aceito no CRAN.

## Pendências e decisões

- **Ordenação PCPS na vignette `whole-tree`.** Pedido em 2026-11-08: além de PCPS 1 contra o
  ambiente, mostrar a ordenação das comunidades no plano PCPS 1 × PCPS 2 com o ambiente como
  tamanho (ou cor) dos círculos, em R base ou em ggplot2 (só na vignette, via `Suggests`).
  Ver se `pcps_test()` guarda os escores do segundo eixo e as posições das espécies (biplot).
- **Imigração na grade.** Decidido em 2026-11-07: argumento `immigration` em
  `simulate_metacommunity()` e `mc_dynamics()`, com `"poisson"` como padrão (código original,
  sementes antigas iguais) e `"multinomial"` conservando os dispersores, como o mcomsimr. No
  cenário da diferença (grade, dispersão 0,001, extirpação 0,02, alfa contando sítios vazios):
  1,51 com Poisson, 1,61 com multinomial e 1,61 no mcomsimr. O multinomial custa o mesmo em
  dispersão baixa e até 20 vezes mais em dispersão alta (10 x 10 sítios, 50 espécies).
- **Erros no artigo de Thompson et al. (2020).** Tabela 1 (beta), faixa de dispersão,
  conversão do ambiente em postos e metade dos imigrantes no cenário de
  competição-colonização não batem com o código. Vale citar isso no artigo do pacote.
- **Dependências do hub.** Decidido em 2026-11-07: `Imports` só com `ape` e pacotes base
  (`graphics`, `grDevices`, `parallel`, `stats`, `utils`); todo o resto em `Suggests`, chamado
  com `requireNamespace()`. O `purrr`, que só dava a barra de progresso, foi para `Suggests`.
- **Código em MATLAB.** Os suplementos não trazem licença. Reimplementar a partir das
  equações e citar; se for copiar a estrutura do código, falar com os autores.
- **Protocolo estático.** O texto de Peres-Neto et al. (2012) dá `h` entre 0,5 e 1 e usa
  **L** como abundância, sem Poisson nos trechos lidos (a etapa 1 não veio). O pcps-sim
  usa `h` entre 0 e 30 e Poisson, e o manuscrito atribui isso a Peres-Neto et al.
  (2012): conferir essa citação, talvez a fonte seja Dray & Legendre (2008).
- **Appendix C.** Não dá a distribuição de `a` nem a escala de `b`; diz que a presença
  vem de um valor "normally generated", o que gera muito ruído (a legenda da Fig. 4
  comenta o ruído). `simulate_logistic()` usa uniforme por padrão e aceita
  `draw = "normal"`. O cenário III é ambíguo. O código de Ives & Helmus (2011,
  `PGLMM_model_simulations.m`, em `~/Downloads/3567789/`) esclarece parte disso: `a` é
  N(0; 0,5), `b` tem desvio 1 na amostra de espécies, a presença vem de `rand` (uniforme),
  e a redução U(0,5; 1) entra uma vez por gradiente, enquanto o Appendix C a aplica uma
  vez. `simulate_logistic()` adotou esses valores na 0.5.1, e a redução fica uma só, como
  no Appendix C. O resultado da Tabela 5 está no artigo de validação.
- **Validar `simulate_logistic()` com PGLMM.** Feito em 2026-11-07, em
  `data-raw/validation_pglmm.R` (phyr 1.1.3, só no script; árvore balanceada de 32 espécies,
  40 sítios, 100 réplicas por cenário, 900 ajustes sem falha). Modelo
  `y ~ x + (1 | sp) + (x | sp__)`, binomial. Sem a redução por sítio, o PGLMM recupera a
  simulação: fração filogenética da variância das inclinações 0,00, 0,48 e 1,00 para lambda
  0, 0,5 e 1, igual à obtida das inclinações verdadeiras; teste da inclinação filogenética com
  5-6% de rejeição em lambda 0, 71% em 0,5 e 100% em 1; variâncias 10-15% abaixo, viés
  conhecido do ajuste binário. Com a redução padrão U(0,5; 1), que o modelo não tem, as
  variâncias caem à metade, mas a fração e as rejeições não mudam; um intercepto aleatório de
  sítio não corrige. `slope_sd` é o desvio amostral das inclinações (padronizadas por
  `scale()`), não o do processo: com lambda 1 nesta árvore, a variância do processo é 1,48.
  Isso ficou explícito na ajuda.
- **treeshapesim.** O artigo de integração com o comdynamite usa `species_traits()`,
  `competition_matrix()`, `landscape_implicit()` e `simulate_metacommunity()`: mudar
  assinatura ou saída delas quebra o site do treeshapesim. Quando `prune_tips()` estiver
  na `main` de lá, pôr o link do artigo no README ou em `evaluating-methods`. O
  comdynamite fica só em `Config/Needs/website` no treeshapesim, sem dependência circular.
- **Sessão anterior.** A sessão que começou o pacote não foi avisada por aqui de que o
  desenvolvimento mudou de lugar. Commits novos dela vão divergir do `main`.
- **pcps-sim.** Se a submissão for duplo-cega, avaliar se o `NEWS.md` pode manter "Code
  from the pcps-sim project". A vignette `evaluating-methods` compara permutações
  livres e do toro, tema próximo do piloto de validade sob H0 daquele projeto.
- **Resultados pré-calculados.** Os `.rds` em `vignettes/` vêm de
  `data-raw/vignette_ecology.R` e `data-raw/vignette_method.R`. Se o modelo mudar,
  é preciso rodar os scripts de novo e revisar os números citados no texto.
- **Páginas órfãs do site.** Resolvido em 2026-11-07: o deploy passou a `clean: true` (com
  `clean-exclude: .nojekyll`), que apaga do `gh-pages` o que o site não gera mais. Na
  comparação com um site construído localmente, sobravam só 5 figuras antigas de artigos.
- **`data-raw/hex_dynamite.R`.** Removido em 2026-11-07; o logo atual (opção F1) vem de
  `data-raw/hex_sticker.R`.

## Referências citadas neste roadmap

- Desdevises, Y., Legendre, P., Azouzi, L., & Morand, S. (2003). *Evolution*, 57, 2647–2652. https://doi.org/10.1111/j.0014-3820.2003.tb01508.x
- Drury, J., Clavel, J., Manceau, M., & Morlon, H. (2016). *Systematic Biology*, 65, 700–710. https://doi.org/10.1093/sysbio/syw020
- Leibold, M. A., & Chase, J. M. (2018). *Metacommunity Ecology*. Monographs in Population Biology 59. Princeton University Press. https://doi.org/10.2307/j.ctt1wf4d24

Conferidas no Crossref em 2026-11-05.

- Benício, R. A., et al. (2021). *Journal of Biogeography*, 48, 2211–2225. https://doi.org/10.1111/jbi.14145
- ter Braak, C. J. F. (2019). *Methods in Ecology and Evolution*, 10, 1962–1971. https://doi.org/10.1111/2041-210X.13278
- Caetano, D. S., et al. (2018). *Evolution*, 72, 2308–2324. https://doi.org/10.1111/evo.13602
- Csilléry, K., et al. (2012). *Methods in Ecology and Evolution*, 3, 475–479. https://doi.org/10.1111/j.2041-210x.2011.00179.x
- Dray, S., & Legendre, P. (2008). *Ecology*, 89, 3400–3412. https://doi.org/10.1890/08-0349.1
- Clavel, J., et al. (2015). *Methods in Ecology and Evolution*, 6, 1311–1319. https://doi.org/10.1111/2041-210x.12420
- Duarte, L. D. S., et al. (2016). *Methods in Ecology and Evolution*, 7, 937–946. https://doi.org/10.1111/2041-210x.12547
- Duarte, L. D. S., et al. (2023). *Ecography*, 2023. https://doi.org/10.1111/ecog.06544
- Goldberg, E. E., et al. (2011). *Systematic Biology*, 60, 451–465. https://doi.org/10.1093/sysbio/syr046
- Goldberg, E. E., & Igić, B. (2012). *Evolution*, 66, 3701–3709. https://doi.org/10.1111/j.1558-5646.2012.01730.x
- Guzman, L. M., et al. (2022). *Ecology*, 103. https://doi.org/10.1002/ecy.3683
- Hubert, N., Calcagno, V., Etienne, R. S., & Mouquet, N. (2015). *Ecology Letters*, 18, 864–881. https://doi.org/10.1111/ele.12458
- Harvey, M. G., & Rabosky, D. L. (2017). *Methods in Ecology and Evolution*, 9, 984–993. https://doi.org/10.1111/2041-210x.12949
- Ives, A. R., & Helmus, M. R. (2011). *Ecological Monographs*, 81, 511–525. https://doi.org/10.1890/10-1264.1
- Legendre, P. (2008). *Journal of Plant Ecology*, 1, 3–8. https://doi.org/10.1093/jpe/rtm001
- Januario, M., Pinsky, M. L., & Rabosky, D. L. (2025). *Ecology Letters*, 28, e70021. https://doi.org/10.1111/ele.70021
- Leibold, M. A., Govaert, L., Loeuille, N., De Meester, L., & Urban, M. C. (2022). *Annual Review of Ecology, Evolution, and Systematics*, 53, 299–326. https://doi.org/10.1146/annurev-ecolsys-102220-024934
- Leibold, M. A., Urban, M. C., De Meester, L., Klausmeier, C. A., & Vanoverbeke, J. (2019). *Proceedings of the National Academy of Sciences*, 116, 2612–2617. https://doi.org/10.1073/pnas.1808615116
- Leibold, M. A., et al. (2010). *Ecology Letters*, 13, 1290–1299. https://doi.org/10.1111/j.1461-0248.2010.01523.x
- Matzke, N. J. (2014). *Systematic Biology*, 63, 951–970. https://doi.org/10.1093/sysbio/syu056
- McPeek, M. A. (2017). *Evolutionary Community Ecology*. Monographs in Population Biology. Princeton University Press. https://doi.org/10.23943/princeton/9780691088778.001.0001
- Morlon, H., et al. (2016). *Methods in Ecology and Evolution*, 7, 589–597. https://doi.org/10.1111/2041-210x.12526
- O'Meara, B. C., Ané, C., Sanderson, M. J., & Wainwright, P. C. (2006). *Evolution*, 60, 922–933. https://doi.org/10.1111/j.0014-3820.2006.tb01171.x
- Nakamura, G., et al. (2023). *Journal of Biogeography*, 51, 533–543. https://doi.org/10.1111/jbi.14610
- Pontarp, M., et al. (2019). *Methods in Ecology and Evolution*, 10, 450–460. https://doi.org/10.1111/2041-210x.13129
- Prangle, D., et al. (2014). *Statistical Applications in Genetics and Molecular Biology*, 13, 67–82. https://doi.org/10.1515/sagmb-2013-0012
- Pudlo, P., et al. (2015). *Bioinformatics*, 32, 859–866. https://doi.org/10.1093/bioinformatics/btv684
- Rabosky, D. L., & Goldberg, E. E. (2015). *Systematic Biology*, 64, 340–355. https://doi.org/10.1093/sysbio/syu131
- Rabosky, D. L., & Goldberg, E. E. (2017). *Evolution*, 71, 1432–1442. https://doi.org/10.1111/evo.13227
- McArdle, B. H., & Anderson, M. J. (2001). *Ecology*, 82, 290–297. https://doi.org/10.1890/0012-9658(2001)082[0290:FMMTCD]2.0.CO;2
- Peres-Neto, P. R., & Kembel, S. W. (2015). *Plant Ecology*, 216, 709–724. https://doi.org/10.1007/s11258-014-0405-0
- Peres-Neto, P. R., et al. (2012). *Ecology*, 93, S14–S30. https://doi.org/10.1890/11-0494.1
- Pigot, A. L., & Etienne, R. S. (2015). *Ecology Letters*, 18, 153–163. https://doi.org/10.1111/ele.12395
- Robert, C. P., et al. (2011). *PNAS*, 108, 15112–15117. https://doi.org/10.1073/pnas.1102900108
- Ruffley, M., et al. (2019). *Ecology and Evolution*, 9, 13218–13230. https://doi.org/10.1002/ece3.5773
- van der Plas, F., et al. (2015). *Ecology*, 96, 1502–1515. https://doi.org/10.1890/14-0454.1
- Urban, M. C., Strauss, S. Y., Pelletier, F., Palkovacs, E. P., Leibold, M. A., Hendry, A. P., De Meester, L., Carlson, S. M., Angert, A. L., & Giery, S. T. (2020). *Proceedings of the National Academy of Sciences*, 117, 17482–17490. https://doi.org/10.1073/pnas.1918960117
- Violle, C., Enquist, B. J., McGill, B. J., Jiang, L., Albert, C. H., Hulshof, C., Jung, V., & Messier, J. (2012). *Trends in Ecology & Evolution*, 27, 244–252. https://doi.org/10.1016/j.tree.2011.11.014
- Vehtari, A., et al. (2017). *Statistics and Computing*, 27, 1413–1432. https://doi.org/10.1007/s11222-016-9696-4
- Wood, S. N. (2010). *Nature*, 466, 1102–1104. https://doi.org/10.1038/nature09319
