# comdynamite: estado e próximos passos

Atualizado em 2026-11-05. Versão 0.1.0 no `DESCRIPTION`, ainda sem tag.
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

## v0.1.0: primeira versão citável

Objetivo: uma versão estável, com DOI, que possa ser citada e instalada.

- [x] Workflow `R-CMD-check.yaml` (macOS, Windows, Linux) e badge no README.
- [x] Medir a cobertura de testes (`covr`): 100% das linhas.
- [x] `inst/CITATION` e `CITATION.cff`. Falta acrescentar o DOI depois do Zenodo.
- [ ] Integração com o Zenodo e DOI de conceito no README.
- [ ] Publicar no r-universe (`diogoprov.r-universe.dev`), como o treeshapesim.
- [x] Versão 0.1.0 no `DESCRIPTION` e no `NEWS.md`.
- [ ] Tag `v0.1.0` e release no GitHub, com o Zenodo já ativado.

**Concluído quando:** o check passar nas três plataformas e o DOI estiver no README.

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
