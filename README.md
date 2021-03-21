# ALGORITMOS-DE-DECODIFICACAO-POR-DECISAO-SUAVE-PARA-CODIGOS-DE-BLOCO-CURTOS
Códigos desenvolvidos para estimar o desempenho dos seguintes algoritmos de decodificação por decisão suave: máxima verossimilhança, Chase II e Pyndiah iterativo para códigos de produto.

Como o intuito do trabalho de conclusão de curso que motivou o desenvolvimento desses códigos MATLAB era avaliar a aplicação de códigos corretores de erros lineares e de bloco a canais com taxas de transmissão limitada a poucos de quilobits por segundo, como o canal acústico submarino, o comprimento dos códigos utilizados foi limitado centenas de bits. Os códigos utilizados foram o código de Hamming, o código binário de Golay, alguns códigos BCH e códigos de produto construidos a partir da combinação dos anteriores.

# 1. Hirarquia de diretórios
#### 1.1 data
Diretório com dados a serem utilizados como input para as funções.
#### 1.2 env
Diretório com variáveis de ambiente para as funções.
#### 1.3 lib
Diretório com funções desenvolvidas para serem utilizadas pelas funções principais.
#### 1.4 results
Diretório com os resultados obtidos. Ao executar uma simulação o resultado é salvo em um arquivo no diretório atual, logo, caso queira guardá-lo de forma organizada, recomenda-se que o mesmo seja movido para o diretório results.

# 2. Parâmetros das funções
As funções principais são as que começam com os prefixos "FER_BER_", "INPUT_" ou "SD_" e os principais parâmetros utilizados são descritos a seguir.
#### 2.1 wethr
número de erros de palavra a serem observados
#### 2.2 output_method
'a' ou 'w' (append ou write)
#### 2.3 radius
valor para raio de Chase
#### 2.4 iterationWeight
vetor com o peso de cada iteração para o algortimo de Pyndiah
#### 2.5 N, K, dmin e t
parâmetros do código a ser simulado
#### 2.6 vEbNo ou "EbNoMin, EbNoMax e step"
vetor com os valores de EbNo a serem simulados ou o início, o fim e o passo para a criação do vetor com os valores de EbNo.

# 3. Gerando gráficos
Os gráficos apresentados no trabalho de conclusão de curso foram gerados a partir do script plotCalls.m que utiliza a função resultsToPlot.m, ambas desenvolvidas pelo autor. Logo, para obter os gráficos basta executar a função plotCalls.m.
