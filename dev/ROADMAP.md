# comdynamite: estado e próximos passos

Atualizado em 2026-11-05. Versão atual: 0.1.0 (tag `v0.1.0`, DOI de conceito `10.5281/zenodo.23993533`).
As versões abaixo são uma sugestão: cada marco tem um critério de conclusão, e a
ordem pode mudar. Tag no formato `vX.Y.Z`, com release no GitHub.

## Feito

| Item | O que entrou |
|---|---|
| Modelo | Thompson et al. (2020) com traits e competição filogenéticos; paisagens implícitas e em grade; permutações do toro; ambiente variável no tempo |
| Montagem por regras | `assemble_community()`: similaridade limitante, filtragem e neutra |
| Traits | `traits_ou()` (um ou vários ótimos) e `paint_regimes()`, no lugar de `traits_bm()` |
| Identidade | hex sticker (opção E) e `data-raw/hex_sticker.R` |
| Site | pkgdown em https://provetelab.org/comdynamite/, com build automático no `main` |
| Vignettes | `comdynamite`, `phylogenetic-structure` e `evaluating-methods` |
| Documentação | `\value` e exemplos em todas as funções exportadas; `Description` atualizada |

## v0.1.0: primeira versão citável (concluída em 2026-11-05)

CI em macOS, Windows e Linux; 100% de cobertura de linhas; `inst/CITATION` e
`CITATION.cff`; release no GitHub com DOI do Zenodo. Falta confirmar o pacote no
r-universe: o registro foi feito, e o build depende do ciclo do serviço.

## v0.2.0: validação do modelo

Objetivo: mostrar que a implementação reproduz o modelo publicado.

- [ ] Comparar com o `mcomsimr` nos mesmos parâmetros, levando em conta as diferenças
      de convenção (`sigma`, `alpha[j, k]`). Comparar distribuições, não réplicas.
- [ ] Reproduzir qualitativamente os padrões de Thompson et al. (2020).
- [ ] Testes de propriedades com resultado analítico (por exemplo, o equilíbrio de
      Beverton-Holt com uma espécie, já usado no exemplo de `mc_dynamics()`).
- [ ] Vignette ou seção "Validation" no site.

**Concluído quando:** houver um script em `data-raw/` que reproduza a comparação.

## v0.3.0: fluxo de estudos de simulação

Objetivo: reduzir o código que o usuário repete. As duas vignettes novas precisaram de
laços e funções escritos à mão, o que indica o que falta no pacote.

- [ ] Função para rodar grades de cenários com réplicas pareadas por árvore e
      paisagem, com semente controlada e paralelismo opcional.
- [ ] Diagnóstico de equilíbrio. Com competição filogenética e dispersão alta, a
      exclusão regional continuou até 3000 passos, e hoje nada avisa o usuário.
- [ ] Saída em formato longo (data frame) a partir de `history`.
- [ ] Benchmark de `mc_dynamics()`. Avaliar código compilado só se o benchmark
      mostrar ganho real, como foi feito no treeshapesim.

**Concluído quando:** as vignettes puderem ser reescritas com essas funções, com bem
menos código.

## v0.4.0: mais processos e integração

Objetivo: ampliar os cenários que o pacote cobre. Escolher pelo uso real.

- [ ] Vários eixos ambientais (hoje há um só).
- [ ] Paisagens irregulares ou lidas de coordenadas, além da grade.
- [ ] Null models restritos para paisagens com bordas. `torus_permutations()` só vale
      em toro, e isso limita o uso com dados reais.
- [ ] Gerador de ambiente temporal fiel ao de Thompson et al. (2020), como alternativa
      a `env_over_time()`.
- [ ] Vignette com árvores do treeshapesim e do gen3sis como pool regional.

**Concluído quando:** os itens escolhidos tiverem testes e exemplo em vignette.

## v0.5.0: geradores estáticos e métodos de metacommunity phylogenetics

Objetivo: reunir os geradores e os testes usados na literatura de validação, para
comparar métodos no mesmo desenho. O código em MATLAB está nos suplementos de
Peres-Neto et al. (2012) e Leibold et al. (2010).

- [ ] Gerador estático de Peres-Neto et al. (2012): resposta gaussiana com traits em
      Brownian motion e, com dois gradientes, produto das respostas ao ambiente e ao
      espaço. Partir de `R/simulate_layerA.R` e `R/simulate_static.R` do pcps-sim, que
      já reproduzem as validações publicadas.
- [ ] Gerador logístico do Appendix C (Ives & Helmus 2011), com os sete cenários de
      [E] e [S].
- [ ] Traduzir de MATLAB para R o framework de Peres-Neto et al. (2012): correlações
      entre traits e ambiente, autovetores ponderados, seleção estocástica e partição
      da variação. Reaproveitar `ade4::fourthcorner` e `adespatial::mem` se servirem.
- [ ] Traduzir `MetacommunityPhylogenetics` de Leibold et al. (2010): regressão
      ponderada por nó, com testes de permutação.
- [ ] Testar a tradução contra o MATLAB em dados fixos, antes de corrigir qualquer
      coisa. Ver os problemas listados em "Pendências e decisões".
- [ ] Paralelismo com o Futureverse e estilo tidyverse, sem pesar nas dependências:
      `future`, `furrr` ou `future.apply` em `Suggests`.

**Concluído quando:** a tradução reproduzir os resultados do MATLAB e os testes
reproduzirem o erro tipo I publicado.

## v0.6.0: mais evolução de traits e modelos nulos

Objetivo: ampliar os modelos de traits e os nulos com dinâmica histórica.

- [ ] Mais opções de evolução de traits, multivariadas e com covariância entre traits,
      tomando o mvMORPH (Clavel et al. 2015) como referência.
- [ ] Modelo nulo dinâmico do DAMOCLES (Pigot & Etienne 2015): colonização e extinção
      local ao longo da árvore. Decidir entre chamar o pacote e reimplementar.

**Concluído quando:** cada opção tiver teste contra o pacote de referência.

## v0.7.0: traits, filogenia e biogeografia

Objetivo: gerar cenários em que história biogeográfica, traits e montagem local estão
ligados, como em Benício et al. (2021), com a verdade conhecida.

- [ ] Fluxo que integra traits, filogenias e biogeografia, seguindo Benício et al.
      (2021).
- [ ] Modelos de Markov com k estados para áreas e traits discretos: DEC e
      variantes do BioGeoBEARS (Matzke 2014), ClaSSE (Goldberg & Igić 2012) e GeoSSE
      (Goldberg et al. 2011).
- [ ] Modelos contínuos: SBEARS, do Herodotools.
- [ ] Integração com o Herodotools (Nakamura et al. 2023): saídas do comdynamite nos
      formatos que ele lê, e as métricas de `calc_tip_based_trait_evo()` como
      resumo. Vignette de exemplo.

**Concluído quando:** um cenário simulado passar pelo Herodotools e recuperar a
história conhecida.

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

- **Escopo.** Os itens de análise da v0.5.0 são métodos, não simulação. Decidir se
  ficam no comdynamite ou em pacote irmão. Tidyverse e Futureverse em `Imports` também
  pesariam: hoje o pacote só importa `ape` e `stats`.
- **Código em MATLAB.** Os suplementos não trazem licença. Reimplementar a partir das
  equações e citar; se for copiar a estrutura do código, falar com os autores.
- **Possíveis erros no MATLAB, a confirmar rodando.** Em `MetacommunityPhylogenetics`,
  nas permutações 2 e 3, `randperm(NumberSpecies)` usa o número de espécies do nó, então
  só sorteia as primeiras colunas da matriz; e `TotalSum` não é recalculado depois
  que `Incidence` muda.
- **Protocolo estático.** `ScriptDemo.m` usa `h` entre 0,5 e 1 e não sorteia contagens;
  o pcps-sim usa `h` entre 0 e 30 e Poisson. Conferir no texto de Peres-Neto et al.
  (2012) antes de definir os padrões da função.
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
- Clavel, J., et al. (2015). *Methods in Ecology and Evolution*, 6, 1311–1319. https://doi.org/10.1111/2041-210x.12420
- Goldberg, E. E., et al. (2011). *Systematic Biology*, 60, 451–465. https://doi.org/10.1093/sysbio/syr046
- Goldberg, E. E., & Igić, B. (2012). *Evolution*, 66, 3701–3709. https://doi.org/10.1111/j.1558-5646.2012.01730.x
- Ives, A. R., & Helmus, M. R. (2011). *Ecological Monographs*, 81, 511–525. https://doi.org/10.1890/10-1264.1
- Leibold, M. A., et al. (2010). *Ecology Letters*, 13, 1290–1299. https://doi.org/10.1111/j.1461-0248.2010.01523.x
- Matzke, N. J. (2014). *Systematic Biology*, 63, 951–970. https://doi.org/10.1093/sysbio/syu056
- Nakamura, G., et al. (2023). *Journal of Biogeography*, 51, 533–543. https://doi.org/10.1111/jbi.14610
- Peres-Neto, P. R., et al. (2012). *Ecology*, 93, S14–S30. https://doi.org/10.1890/11-0494.1
- Pigot, A. L., & Etienne, R. S. (2015). *Ecology Letters*, 18, 153–163. https://doi.org/10.1111/ele.12395
