//Além disso, ao final de cada bloco (após a matriz inversa) deve constar:
//O tempo de execução da triangularização;
//O tempo de cálculo de Y e X nas fases Ly=b e Ux=y
//O valor da norma L2 do resíduo de A.A' = I, isto é, a norma L2 do resíduo de cada coluna da matriz 
//identidade com a coluna correspondente da matriz A.A', onde A'  é o resultado do cálculo de A-1 pelo Método da Fatoração LU.

//Desenvolvolvedores:
	//Anna Caroline Bozzi - acb17 -  GRR20173532
	//João Vitor Gualberto Figueiredo - jvgf17 - GRR20177037

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "matriz.h"
#include "utils.h"
#include <ctype.h>
#include <unistd.h>
#include <stdbool.h>
#include <getopt.h>


int main (int argc, char *argv[])
{

	int tamanho;
	real_t *resi;
	double *tTotal;


	int opcao;
	char *pathfile = NULL;
	int xflag = 0;
	int file_write = 0;
    FILE *file;


	Matriz_t *SL;

	
    while((opcao = getopt(argc, argv, ":po:"))!=-1)
    {
		switch(opcao)
        {
			case 'p': //se passa a opção -p é para fazer fatoração com pivoteamento
				xflag = 1;
                break;

    		case 'o':
    			file_write = 1;
    			pathfile = optarg;
                

                file = fopen(pathfile, "a");
                if (file == NULL)
                {
                    fprintf(stderr, "Erro ao abrir o arquivo.\n");
                    exit (-1);
                }

                break;
		}
	}
    


    char arq = 'x';
	while (arq!=EOF){
		if(arq != ' '){
			SL = lerMatriz();
			
			tamanho = SL->n;
			Matriz_t *original;
			original = alocaMatriz(tamanho);
			copiaMatriz(SL,original);

			resi = (real_t *) malloc(tamanho*sizeof(real_t));

			tTotal = (double *) malloc(tamanho*sizeof(double));

			for(int i=0; i<tamanho;i++){
				resi[i]=0;
				tTotal[i] = 0;
			}

			if (xflag == 1){
				double time = timestamp();
				decomposicaoLUpivoteamento(SL,tTotal);
    			double elapsed = timestamp() - time;
    			tTotal[2] += elapsed;
			}else{
				double time = timestamp();
				decomposicaoLU(SL,tTotal);
    			double elapsed = timestamp() - time;
    			tTotal[2] += elapsed;
			}

			for(int i=0; i<tamanho;i++){
				for(int j=0;j<tamanho;j++){
					if (isinf(SL->A[i][j]) || isnan(SL->A[i][j])){
						printf("Ocorreu operações com 0, o programa foi encerrado\n"); //trocar para fprint (stderr, "Ocorreu operações com 0, o programa foi encerrado\n")
						exit(0);
					}
				}
			}


			normaL2residuo(original,SL,resi);


            //se for uma escrita num arquivo, passa o arquivo já aberto
            if (file_write == 1)
            {
                fprintf(file, "%d\n", tamanho);
			    prnMatriz(original, file_write, file); // original é a original
			    fprintf(file, "#\n");
			    prnMatriz(SL, file_write, file); //SL agora é a inversa
			    fprintf(file, "###########\n");
			    fprintf(file, "# Tempo Triangularização: %.5g milissegundos. \n",tTotal[2]);
			    fprintf(file, "# Tempo cálculo de Y: %.5g milissegundos. \n",tTotal[0]);
			    fprintf(file, "# Tempo cálculo de X: %.5g milissegundos. \n", tTotal[1]);
			    fprintf(file, "# Norma L2 do Residuo:");
			    for(int i=0; i<tamanho;i++)
				    fprintf(file, "%g ",resi[i]);
                fprintf(file, "\n");
			    fprintf(file, "###########\n");
			
            }
            //se esse argumento não for passado, imprima na saída padrão
            else
            {
                fprintf(stdout, "%d\n", tamanho);
			    prnMatriz(original, file_write, file); // original é a original
			    fprintf(stdout, "#\n");
			    prnMatriz(SL, file_write, file); //SL agora é a inversa
			    fprintf(stdout, "###########\n");
			    fprintf(stdout, "# Tempo Triangularização: %.5g milissegundos. \n",tTotal[2]);
			    fprintf(stdout, "# Tempo cálculo de Y: %.5g milissegundos. \n",tTotal[0]);
			    fprintf(stdout, "# Tempo cálculo de X: %.5g milissegundos. \n", tTotal[1]);
			    fprintf(stdout, "# Norma L2 do Residuo:");
			    for(int i=0; i<tamanho;i++)
				    fprintf(stdout, "%g ",resi[i]);
                fprintf(stdout, "\n");
			    fprintf(stdout, "###########\n");
			}
		}
		arq = getchar();
	}
}