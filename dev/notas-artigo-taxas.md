# Notas para o artigo sobre taxas por ramo

Leituras para o artigo (vignette) sobre taxas de evolução de traits que variam entre ramos e o
papel da montagem de metacomunidades. Distingue o que cada fonte afirma da nossa interpretação.
Citações literais conferidas contra o PDF ou o texto extraído. Referências completas na seção
"Referências" do `ROADMAP.md`.

## Pergunta do artigo

A montagem da metacomunidade molda a evolução dos traits e deixa rastro nas taxas estimadas por
ramo? `trait_rates()` e `trait_mk()` dão a verdade conhecida; a validação em
`data-raw/validation_rates.R` e `data-raw/validation_bamm.R` mostra quanto os estimadores
recuperam. Limite central: com uma realização por ramo, a taxa de cada ramo é pouco
identificável (correlação entre mudança verdadeira e taxa verdadeira de 0,29 a 0,40).

## Hubert et al. (2015)

- Revisão dos modelos de especiação em metacomunidades (neutros, de nicho, de trade-offs de
  história de vida).
- Especiação em duas etapas: surgimento contínuo de espécies incipientes e persistência delas por
  um período de transição; estrutura espacial e dispersão regulam as duas.

## Januario, Pinsky & Rabosky (2025)

- Arcabouço para testar se a persistência metapopulacional (ocupação de equilíbrio e
  resiliência) prediz a taxa de especiação nas pontas.
- Peixes demersais: um índice de persistência tem sinal filogenético, mas nenhum prediz a taxa
  de especiação; os autores falam em desacoplamento entre a dinâmica de décadas e a de milhões
  de anos.
- Os dois trabalhos tratam de taxas de diversificação; aqui as taxas são de evolução de traits.

## Leibold, Govaert, Loeuille, De Meester & Urban (2022)

- A evolução adaptativa altera a montagem local e regional, e o efeito depende da
  heterogeneidade ambiental e da sua estrutura espacial. Resultado típico: efeitos de prioridade
  mediados pela evolução (colonizadores precoces se adaptam e alteram o estabelecimento dos que
  chegam depois; monopolização da comunidade).
- Três taxas relativas decidem o desfecho (Fig. 1): colonização, regulação da comunidade e
  evolução. Evolução lenta: quadro ecológico clássico. Evolução rápida em relação à colonização:
  monopolização e, no extremo, radiação adaptativa. Colonização alta: fluxo gênico afoga a
  adaptação local (gene swamping).
- Competição: deslocamento de caráter, convergência que produz coocorrência quase neutra quando
  há mais espécies que nichos, ciclos eco-evolutivos. Interações tróficas: perda de defesas na
  ausência do predador (genetic subsidy).
- Lacuna apontada pelos autores: os modelos macroecológicos não incorporam adequadamente quanto
  a adaptação local e a montagem alteram as taxas de especiação e extinção e, por elas, as
  filogenias e a biota regional.
- Cenários sugeridos: frequência de extinção de manchas (distúrbios reiniciam a montagem; em
  excesso impedem a adaptação), autocorrelação espacial do ambiente, ordem de chegada,
  adaptações secundárias que revertem a monopolização.
- Métodos de partição de variância de traits: estatísticas T de Violle et al. (2012), Lepš et
  al. (2011). As estatísticas T estão em `cati::Tstats()` (CRAN 0.99.6, GPL >= 2; pode ir para
  `Suggests`): T_IP.IC (variância dentro das populações sobre a variância dentro das
  comunidades), T_IC.IR (comunidades sobre o pool regional) e T_PC.PR (entre comunidades sobre o
  pool regional), cada uma com seu modelo nulo (local, regional.ind, regional.pop). Exigem traits
  de indivíduos, com fatores de parcela e de espécie.
- Interpretação nossa: eles falam de taxas de evolução relativas à colonização em escala
  ecológica, não de taxas por ramo numa filogenia.

## McPeek (2017), *Evolutionary Community Ecology*

Leitura feita por subagente; três citações conferidas no texto extraído.

- Mecanismo: a mudança do traço médio é V_z vezes a soma, sobre as espécies que interagem, da
  abundância de cada uma vezes o gradiente de seleção (eq. 3.7). V_z é constante e pequeno; toda
  variação de velocidade vem do contexto da comunidade (quem está presente, abundâncias, traits).
- Não usa o vocabulário de métodos comparativos (BM, sigma², early burst, taxas por ramo). A
  ligação com modelos de taxas por ramo é interpretação nossa: o livro dá o mecanismo, não a
  estatística macroevolutiva.
- Evolução acelera: invasão de uma comunidade (mudança rápida que desacelera até o equilíbrio,
  Fig. 3.2–3.3), mudanças de hábitat ou hospedeiro ("must occur very rapidly, probably over, at
  most, tens of generations", p. 149), peak shifts na diferenciação, ciclos de traits,
  reembaralhamento das comunidades no Quaternário. "Large and changing types of selection should
  be expected most when species are first adapting to a new selection regime" (p. 129).
- Evolução desacelera: equilíbrios evolutivos (estase com seleção direcional nos componentes),
  espécies neutras (até 12 *Enallagma* por cerca de 10 milhões de anos), adição de predador
  (enfraquece a seleção sobre o recurso), fluxo gênico.
- Deslocamento de caráter não é inevitável: exige abundâncias altas, alta produtividade,
  interações fortes; propaga-se de recursos para consumidores e predadores.
- Metacomunidade (cap. 6): dispersão homogeneiza abundâncias e traits; acima de um limiar de
  diferença entre ótimos locais há transição abrupta para baixa dispersão e diferenciação local.
  Os modelos do cap. 6 têm uma espécie em duas ou J manchas, não uma metacomunidade
  multiespecífica com traits evoluindo.
- "By reconstructing trait evolution across multiple interacting clades, hypotheses about
  coevolutionary responses among species can also be tested" (p. 274).

## Cenários de simulação derivados (interpretação nossa)

1. Explosão na colonização: mudança de taxa no ramo em que a linhagem coloniza um novo tipo de
   mancha, decaindo ao longo do ramo; a linhagem residente responde com atraso.
2. Mudança de hábitat como gatilho: classe de taxa rápida em `trait_mk()` (ou mudança em
   `trait_rates()`) nos ramos com transição de hábitat, seguida de estase.
3. Taxa como função do ambiente: taxa por linhagem crescendo com a produtividade ou a força das
   interações nas manchas ocupadas.
4. Posição trófica: recursos mais rápidos que consumidores, consumidores mais rápidos que
   predadores; queda da taxa do recurso quando um predador entra.
5. Estase pontuada por mudanças na composição: taxa do ramo proporcional ao número de entradas e
   saídas de parceiros de interação ao longo dele, contra um BM homogêneo como nulo. É a previsão
   mais defensável de que a montagem deixa rastro nas taxas.
6. Taxas acopladas entre clados que interagem (recurso e consumidor).
7. Dispersão: linhagens generalistas de alta dispersão com taxas baixas; especialistas de baixa
   dispersão divergem; variar a diferença entre ótimos para achar o limiar.
8. Mudanças de taxa por Poisson com intensidade variável no tempo, concentradas em eventos de
   reembaralhamento das comunidades e compartilhadas entre linhagens simpátricas, contra mudanças
   independentes entre linhagens.
9. Frequência de extinção de manchas e autocorrelação espacial do ambiente (Leibold et al. 2022).
10. Variação intraespecífica e adaptação local: gerar indivíduos a partir das abundâncias
    simuladas, com traits em torno da média da espécie ou ajustados ao ambiente do sítio, e
    calcular as estatísticas T com `cati::Tstats()`. Adaptação local forte deve baixar T_IP.IC
    e mudar T_PC.PR em relação ao cenário sem adaptação.

Para os cenários 1, 2, 5 e 8 o pacote já tem quase tudo: `simulate_metacommunity()` dá a
história de ocupação e `trait_rates(shift_edges = )` ou `trait_mk(edge_rates = )` aceitam taxas
definidas por ramo. Falta a função que traduz a história da metacomunidade em taxas por ramo.
