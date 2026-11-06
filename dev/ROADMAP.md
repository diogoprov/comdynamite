# comdynamite: estado e próximos passos

Atualizado em 2026-11-05. Versão atual: 0.1.0 (tag `v0.1.0`, DOI de conceito `10.5281/zenodo.23993533`).
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
e um prefixo a definir para analisar). O título e o `Description` mudam quando os
primeiros métodos entrarem.

## Feito

| Item | O que entrou |
|---|---|
| Modelo | Thompson et al. (2020) com traits e competição filogenéticos; paisagens implícitas e em grade; permutações do toro; ambiente variável no tempo |
| Montagem por regras | `assemble_community()`: similaridade limitante, filtragem e neutra |
| Traits | `traits_ou()` (um ou vários ótimos) e `paint_regimes()`, no lugar de `traits_bm()` |
| Identidade | hex sticker (opção E) e `data-raw/hex_sticker.R` |
| Site | pkgdown em https://provetelab.org/comdynamite/, com build automático no `main` |
| Simulador estático | `simulate_static()`, na versão de desenvolvimento |
| Vignettes | `comdynamite`, `phylogenetic-structure` e `evaluating-methods` |
| Documentação | `\value` e exemplos em todas as funções exportadas; `Description` atualizada |

## v0.1.0: primeira versão citável (concluída em 2026-11-05)

CI em macOS, Windows e Linux; 100% de cobertura de linhas; `inst/CITATION` e
`CITATION.cff`; release no GitHub com DOI do Zenodo.

- [x] Submissão ao r-universe: `comdynamite` incluído em `packages.json`, no repositório
      `diogoprov/diogoprov.r-universe.dev`.
- [ ] Confirmar o build em https://diogoprov.r-universe.dev/comdynamite. Depois, pôr o
      badge e a instalação pelo r-universe no README, como no treeshapesim.

## v0.2.0: validação e geradores estáticos

Objetivo: mostrar que a implementação reproduz o modelo publicado e ter uma linha de
base sem dinâmica.

- [x] `simulate_static()`: resposta gaussiana com contagens de Poisson (Peres-Neto
      et al. 2012).
- [x] Dois gradientes em `simulate_static()`: produto das respostas ao ambiente e ao
      espaço, como em `ScriptDemo.m`.
- [x] Gerador logístico do Appendix C (Ives & Helmus 2011), com os sete cenários de
      [E] e [S]: `simulate_logistic()` e `logistic_gradients()`.
- [ ] Comparar com o `mcomsimr` nos mesmos parâmetros, levando em conta as diferenças
      de convenção (`sigma`, `alpha[j, k]`). Comparar distribuições, não réplicas.
- [ ] Reproduzir qualitativamente os padrões de Thompson et al. (2020).
- [ ] Testes de propriedades com resultado analítico (por exemplo, o equilíbrio de
      Beverton-Holt com uma espécie, já usado no exemplo de `mc_dynamics()`).
- [ ] Vignette ou seção "Validation" no site.

**Concluído quando:** houver um script em `data-raw/` que reproduza a comparação.

## v0.3.0: fluxo de estudos de simulação

Objetivo: reduzir o código que o usuário repete. As vignettes precisaram de laços e
funções escritos à mão, o que indica o que falta no pacote.

- [ ] Função para rodar grades de cenários com réplicas pareadas por árvore e
      paisagem e semente controlada, com paralelismo pelo Futureverse (`future`,
      `furrr` ou `future.apply` em `Suggests`).
- [ ] Diagnóstico de equilíbrio. Com competição filogenética e dispersão alta, a
      exclusão regional continuou até 3000 passos, e hoje nada avisa o usuário.
- [ ] Saída em formato longo (tibble) a partir de `history`.
- [ ] Benchmark de `mc_dynamics()`. Avaliar código compilado só se o benchmark
      mostrar ganho real, como foi feito no treeshapesim. O ABC da v0.5.0 depende disso.
- [ ] `summary_stats()`: as mesmas estatísticas para dados simulados e observados
      (alfa, beta, gama, ocupação, SAD, SES MPD e MNTD, decaimento com a distância,
      relação entre traits e ambiente).
- [ ] Piloto de identificabilidade, com poucos parâmetros: quanto dá para recuperar a
      partir dessas estatísticas? É o ponto de decisão para a v0.5.0, e o desenho deve
      ser combinado antes de rodar.

**Concluído quando:** as vignettes puderem ser reescritas com essas funções, com bem
menos código.

## v0.4.0: base do hub

Objetivo: uma interface comum, para que simulação e métodos troquem dados entre si e
com outros pacotes. As versões seguintes dependem desta.

- [ ] Estrutura de dados comum: comunidade, árvore, traits, ambiente e coordenadas,
      com checagem de nomes. Serve para dados simulados e reais.
- [ ] Árvore e distribuições como entrada: pool regional por sítio e estado inicial
      vindos de fora. Cobre gen3sis, modelos de SSE e reconstruções biogeográficas.
- [ ] `species_traits()` aceitar traits vindos de fora, por exemplo do mvMORPH.
- [ ] Registrar a história dos traits ao longo dos ramos, se hoje só houver as pontas.
- [ ] Integração com o Herodotools (Nakamura et al. 2023): saídas nos formatos que ele
      lê e as métricas de `calc_tip_based_trait_evo()` como resumo, com vignette.
      Com a história conhecida, a simulação pode validar esses métodos.

**Concluído quando:** um cenário simulado passar pelo Herodotools e recuperar a
história conhecida.

## v0.5.0: inferência por ABC

Objetivo: estimar parâmetros e comparar modelos a partir de dados observados. Depende do
piloto da v0.3.0 e da estrutura comum da v0.4.0. Exemplos de uso de ABC em montagem de
comunidades: Pontarp et al. (2019).

- [ ] Tabela de referência: sortear parâmetros dos priors e simular, com a árvore, o
      ambiente e as coordenadas observados e com detecção imperfeita.
- [ ] Estimação e escolha de modelos com o pacote `abc` (Csilléry et al. 2012).
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

**Concluído quando:** uma vignette levar dados reais do ajuste à avaliação, com
cobertura medida em dados simulados.

## v0.6.0: métodos de metacommunity phylogenetics

Objetivo: levar para o R os métodos que só existem em MATLAB. O código está nos
suplementos de Peres-Neto et al. (2012) e Leibold et al. (2010), e nenhum dos dois
métodos está em pacote de R.

- [ ] Peres-Neto et al. (2012): correlações entre traits e ambiente, autovetores
      ponderados, seleção estocástica e partição da variação. Reaproveitar
      `ade4::fourthcorner` e `adespatial::mem` se servirem.
- [ ] Leibold et al. (2010), `MetacommunityPhylogenetics`: regressão ponderada por nó,
      com testes de permutação.
- [ ] Estilo tidyverse na entrada e na saída, com álgebra de matrizes por dentro;
      permutações em paralelo com o Futureverse.
- [ ] Testar a tradução contra o MATLAB em dados fixos, antes de corrigir qualquer
      coisa. Ver os possíveis erros em "Pendências e decisões".
- [ ] Vignette com os métodos aplicados a dados simulados e reais.

**Concluído quando:** a tradução reproduzir os resultados do MATLAB e os testes
reproduzirem o erro tipo I publicado.

## v0.7.0: evolução de traits e modelos nulos

Objetivo: ampliar os modelos de traits e os nulos com dinâmica histórica, chamando os
pacotes de referência em vez de reimplementá-los.

- [ ] Mais opções de evolução de traits, multivariadas e com covariância entre traits,
      pelo mvMORPH (Clavel et al. 2015).
- [ ] Modelo nulo dinâmico do DAMOCLES (Pigot & Etienne 2015): colonização e extinção
      local ao longo da árvore.

**Concluído quando:** cada opção tiver teste e exemplo em vignette.

## v0.8.0: traits, filogenia e biogeografia

Objetivo: transformar o desenho de Benício et al. (2021) em um fluxo do pacote e testar,
com a verdade conhecida, se ele separa tempo para especiação, taxas de diversificação
e conservação de nicho. Scripts e dados do artigo estão no Dryad
(https://doi.org/10.5061/dryad.sf7m0cg65).

Simular o cenário:

- [ ] Árvore e regiões geradas juntas, por exemplo com um simulador de GeoSSE, em cenários
      que variem especiação, extinção, dispersão e tempo de colonização por região.
- [ ] Traits discretos e filtragem por eles, como os modos reprodutivos. Hoje o modelo
      só tem traits contínuos.
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
- [ ] Null models restritos para paisagens com bordas. `torus_permutations()` só vale
      em toro, e isso limita o uso com dados reais.
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

- **Dependências do hub.** Hoje o pacote só importa `ape` e `stats`. Tidyverse, Futureverse
  e pacotes de biogeografia pesam em `Imports`; preferir `Suggests` sempre que der.
- **Código em MATLAB.** Os suplementos não trazem licença. Reimplementar a partir das
  equações e citar; se for copiar a estrutura do código, falar com os autores.
- **Possíveis erros no MATLAB, a confirmar rodando.** Em `MetacommunityPhylogenetics`,
  nas permutações 2 e 3, `randperm(NumberSpecies)` usa o número de espécies do nó, então
  só sorteia as primeiras colunas da matriz; e `TotalSum` não é recalculado depois
  que `Incidence` muda.
- **Protocolo estático.** O texto de Peres-Neto et al. (2012) dá `h` entre 0,5 e 1 e usa
  **L** como abundância, sem Poisson nos trechos lidos (a etapa 1 não veio). O pcps-sim
  usa `h` entre 0 e 30 e Poisson, e o manuscrito atribui isso a Peres-Neto et al.
  (2012): conferir essa citação, talvez a fonte seja Dray & Legendre (2008).
- **Appendix C.** Não dá a distribuição de `a` nem a escala de `b`; diz que a presença
  vem de um valor "normally generated", o que gera muito ruído (a legenda da Fig. 4
  comenta o ruído). `simulate_logistic()` usa uniforme por padrão e aceita
  `draw = "normal"`. O cenário III é ambíguo. Conferir em Ives & Helmus (2011).
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
- **`data-raw/hex_dynamite.R`.** Sobrescreve o `logo.png` oficial se for rodado. Remover
  ou mudar o arquivo de saída.

## Referências citadas neste roadmap

Conferidas no Crossref em 2026-11-05.

- Benício, R. A., et al. (2021). *Journal of Biogeography*, 48, 2211–2225. https://doi.org/10.1111/jbi.14145
- Caetano, D. S., et al. (2018). *Evolution*, 72, 2308–2324. https://doi.org/10.1111/evo.13602
- Csilléry, K., et al. (2012). *Methods in Ecology and Evolution*, 3, 475–479. https://doi.org/10.1111/j.2041-210x.2011.00179.x
- Clavel, J., et al. (2015). *Methods in Ecology and Evolution*, 6, 1311–1319. https://doi.org/10.1111/2041-210x.12420
- Goldberg, E. E., et al. (2011). *Systematic Biology*, 60, 451–465. https://doi.org/10.1093/sysbio/syr046
- Goldberg, E. E., & Igić, B. (2012). *Evolution*, 66, 3701–3709. https://doi.org/10.1111/j.1558-5646.2012.01730.x
- Harvey, M. G., & Rabosky, D. L. (2017). *Methods in Ecology and Evolution*, 9, 984–993. https://doi.org/10.1111/2041-210x.12949
- Ives, A. R., & Helmus, M. R. (2011). *Ecological Monographs*, 81, 511–525. https://doi.org/10.1890/10-1264.1
- Leibold, M. A., et al. (2010). *Ecology Letters*, 13, 1290–1299. https://doi.org/10.1111/j.1461-0248.2010.01523.x
- Matzke, N. J. (2014). *Systematic Biology*, 63, 951–970. https://doi.org/10.1093/sysbio/syu056
- Morlon, H., et al. (2016). *Methods in Ecology and Evolution*, 7, 589–597. https://doi.org/10.1111/2041-210x.12526
- Nakamura, G., et al. (2023). *Journal of Biogeography*, 51, 533–543. https://doi.org/10.1111/jbi.14610
- Pontarp, M., et al. (2019). *Methods in Ecology and Evolution*, 10, 450–460. https://doi.org/10.1111/2041-210x.13129
- Rabosky, D. L., & Goldberg, E. E. (2015). *Systematic Biology*, 64, 340–355. https://doi.org/10.1093/sysbio/syu131
- Rabosky, D. L., & Goldberg, E. E. (2017). *Evolution*, 71, 1432–1442. https://doi.org/10.1111/evo.13227
- Peres-Neto, P. R., et al. (2012). *Ecology*, 93, S14–S30. https://doi.org/10.1890/11-0494.1
- Pigot, A. L., & Etienne, R. S. (2015). *Ecology Letters*, 18, 153–163. https://doi.org/10.1111/ele.12395
- Vehtari, A., et al. (2017). *Statistics and Computing*, 27, 1413–1432. https://doi.org/10.1007/s11222-016-9696-4
- Wood, S. N. (2010). *Nature*, 466, 1102–1104. https://doi.org/10.1038/nature09319
