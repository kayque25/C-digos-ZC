# 🏖️ Análise da Zona Costeira Paraense — Setor 5

Este repositório contém os scripts utilizados no projeto de **extração e análise da linha de costa** no Setor 5 do estado do Pará. A análise foi conduzida com base em dados do Google Earth Engine (GEE), métricas do DSAS (Digital Shoreline Analysis System) e cálculos de indicadores costeiros como **NSM**, **EPR**, **LRR**, bem como estimativas de erro e incerteza associadas.

## 📂 Estrutura do Projeto

project_ZC/<br>
├── Código_GEE.txt # Script GEE para extração de linhas de costa <br>
├── 2_CALCULO_ERROS_RMSE.py # Cálculo do erro RMSE das linhas extraídas<br>
├── 3_ERROS_INCERTEZAS_EPR_LRR.py # Incertezas dos indicadores EPR e LRR<br>
├── 4_Compilação_DSAS.py # Compilação dos resultados do DSAS<br>
├── 5_SETOR_GERAL_DSAS_NSM_EPR_LRR.py # Processamento completo dos indicadores para o setor 5<br>
├── 6_MATRIZ.py # Geração da matriz final com resultados<br>
├── 7_Dados_Matriz_DSAS.py # Organização final dos dados para análise


## 🌍 Objetivo

Monitorar e analisar a dinâmica costeira do litoral paraense (Setor 5), por meio de:

- Extração multitemporal de linhas de costa via GEE.
- Processamento no ArcGIS com o DSAS.
- Cálculo dos indicadores de mudança costeira.
- Avaliação de erros e incertezas associadas.
- Organização em matrizes e relatórios para análises espaciais.

## ⚙️ Tecnologias Utilizadas

- **Python 3** — Processamento dos dados e análises quantitativas.
- **Google Earth Engine (GEE)** — Extração das linhas de costa via script JavaScript.
- **DSAS (Digital Shoreline Analysis System)** — Plugin do ArcGIS para análise de mudanças costeiras.
- **ArcGIS / QGIS** — Suporte ao pré-processamento espacial.

## 📊 Indicadores Costeiros

- **NSM (Net Shoreline Movement)** — Deslocamento líquido da linha de costa.
- **EPR (End Point Rate)** — Taxa de variação entre primeira e última linha.
- **LRR (Linear Regression Rate)** — Taxa baseada em regressão linear.
- **RMSE** — Erro quadrático médio para validação de linhas.
- **Incertezas** — Propagação de erros nos cálculos de EPR e LRR.

## 🚀 Como Executar

1. Clone este repositório:
   ```bash
   git clone https://github.com/seu-usuario/nome-do-repo.git
   
2. Instale as dependências:
   ```bash
   pip install pandas numpy
  
3. Execute os scripts na ordem recomendada:
   ```bash
   Código_GEE.txt (no GEE)
   2_CALCULO_ERROS_RMSE.py
   3_ERROS_INCERTEZAS_EPR_LRR.py
   4_Compilação_DSAS.py
   5_SETOR_GERAL_DSAS_NSM_EPR_LRR.py
   6_MATRIZ.py
   7_Dados_Matriz_DSAS.py

📎 Créditos
Este projeto foi desenvolvido por Kayque da Silva Dias, geógrafo e mestrando na UFPA, com experiência em sensoriamento remoto e análise costeira na região amazônica e Rafael Alexandre Alves Menezes, geógrafo, com atuação voltada à análise ambiental e estudos costeiros na Amazônia.

## 👥 Autores

- **Kayque da Silva Dias**, geógrafo e mestrando na UFPA, com experiência em sensoriamento remoto, geoprocessamento e análise costeira na região amazônica.  
  [Currículo Lattes](http://lattes.cnpq.br/5961292748412062)

- **Rafael Alexandre Alves Menezes**, geógrafo, com atuação voltada à análise ambiental e estudos costeiros na Amazônia.  
  [Currículo Lattes](http://lattes.cnpq.br/4681123810496065)



---

