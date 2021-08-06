//calcular a matriz inversa de várias matrizes quadradas utilizando o Método de Fatoração LU OK
//Além disso, ao final de cada bloco (após a matriz inversa) deve constar:
//O tempo de execução da triangularização;
//O tempo de cálculo de Y e X nas fases Ly=b e Ux=y
//O valor da norma L2 do resíduo de A.A' = I, isto é, a norma L2 do resíduo de cada coluna da matriz 
//identidade com a coluna correspondente da matriz A.A', onde A'  é o resultado do cálculo de A-1 pelo Método da Fatoração LU.

//Desenvolvolvedores:
	//Anna Caroline Bozzi - acb17 -  GRR20173532
	//João Vitor Gualberto Figueiredo - jvgf17 - GRR20177037

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "matriz.h"
#include <float.h>
#include "utils.h"
#include <stdbool.h>

// Exibe a matriz na saída padrão
void prnMatriz (Matriz_t *matriz, int file_write, FILE *file)
{
    if (file_write == 1){
        for (int lin = 0; lin < matriz->n; lin++){
    	    for (int col=0; col < matriz->n ; col++){
      		    fprintf(file, "%f ", matriz->A[lin][col]);
   		    }
    	    fprintf(file, "\n");
    	}
    }else{
        for (int lin = 0; lin < matriz->n; lin++){
    	    for (int col=0; col < matriz->n ; col++){
      		    fprintf(stdout, "%f ", matriz->A[lin][col]);
   		    }
    	    fprintf(stdout, "\n");
  		}
  	}
}


//Leitura de matrizes da Entrada padrão(stdin)
//retorna uma matriz.
//caso erro retorna NULL
Matriz_t *lerMatriz ()
{

 	int n;

	scanf("%d", &n);
 
 	//aloca a matriz, e atribui a dimensão a estrutura de dados escolhida
  	Matriz_t *matriz = alocaMatriz (n);
  	matriz->n = n;  //dimensão

  	//leitura dos elementos e atribuição
  	for (int lin = 0; lin<n; lin++){
    	for (int col=0; col<n ; col++){
      	scanf("%f", &matriz->A[lin][col]);
		  }
  	}

	//se leu e alocou corretamente retorna a matriz  
  	if (matriz){
    	return matriz;
    }
  	//se não leu corretamente retorna nulo
    return NULL;
}

void copiaMatriz(Matriz_t *SL, Matriz_t *original){

	original->n = SL->n;

	for(int linha = 0; linha< SL->n;linha++){
		for(int coluna = 0; coluna < SL->n; coluna ++){
			original->A[linha][coluna] = SL->A[linha][coluna];
		}
	}
}

//alocação de matriz
Matriz_t* alocaMatriz (unsigned int n){

	Matriz_t *matriz = (Matriz_t *) malloc(sizeof(Matriz_t));
  	if ( matriz ) 
  	{
    	matriz->A = (real_t**) malloc(n * sizeof(real_t*)); //aloca coluna
    	for(int i=0; i<n; i++){
      		matriz->A[i] = (real_t *) malloc(n * sizeof(real_t));
    	}
  	}
	return (matriz);

}


//####################INÍCIO DAS OPERAÇÕES COM MATRIZES####################

//##CALCULA O DETERMINANTE PARA AVALIAR SE A MATRIZ POSSUI INVERSA OU NÃO##
//REMOÇÃO DO SEGUNDO FOR E DO IF
int calculaDeterminante(Matriz_t *matriz, int determ){

	for(int i=0;i<matriz->n;i++){			
		determ = matriz->A[i][i]*determ;
	}

	if(determ == 0) //se determinante 0 a matriz não possui inversa
		return 0;
	else
		return 1;
}

//####################CALCULA A NORMA DO RESIDUO####################
real_t *normaL2residuo(Matriz_t *original, Matriz_t *inversa, real_t *resi){
	//multiplicar a original pela inversa ai obtem uma matriz
	
	real_t aux;
	Matriz_t* obtida;
	bool res = true;
	int tam = original->n;
	obtida = alocaMatriz(tam);
	
	double soma = 0;


	for (int linha = 0;linha<tam;linha++){
		for(int coluna = 0;coluna<tam;coluna++){
			for(int i = 0;i<tam;i++){
				aux = aux + (original->A[linha][i]*inversa->A[i][coluna]);
			}
			obtida->A[linha][coluna] = aux;
			aux = 0;
		}
	}

	Matriz_t *ident = alocaMatriz(tam);

  	for (int linha=0; linha<tam; linha++){
    	for (int coluna = 0; coluna<tam; coluna++){
      		if(linha == coluna){
        		ident->A[linha][coluna] = 1;
      		}else{
        		ident->A[linha][coluna] = 0;
      		}
    	}
  	}

	for (int linha = 0;linha<tam;linha++){
		for(int coluna = 0;coluna<tam;coluna++){
			if((obtida->A[linha][coluna]) != (ident->A[linha][coluna])){
				res = false;
			}
		}
	}

	if (res==false){
		for (int i = 0; i<tam;i++){
			soma = 0;
			for(int j =0 ;j<tam;j++){
				soma += fabs(obtida->A[j][i] - ident->A[j][i]);
			}
			soma = sqrt(soma);
			resi[i] = soma;
		}
	}else

		return 0;
}


//####################SOLUÇÃO SISTEMA TRIANGULAR INFERIOR(LY=B)####################
real_t* triangInferior(Matriz_t *L, real_t *x, double *ty){
 
  int n=L->n; //tamanho da matriz L -> decomposição Lu

  real_t *y = (real_t *) malloc(n * sizeof(real_t)); // vetor de solução Y
  real_t soma;
  
  y[0] = x[0]/L->A[0][0];	//x é a coluna da matriz identidade correspondente a solução
  for(int linha = 1; linha< n; linha++){
    soma = 0;
    for(int coluna = 0; coluna <=linha;coluna++){
      soma += L->A[linha][coluna]*y[coluna];
    }
    y[linha] = (x[linha]-soma)/(L->A[linha][linha]);
  }


  return y;
} //OK


//####################SOLUÇÃO SISTEMA TRIANGULAR SUPERIOR(UX=Y)####################
real_t* triangSuperior(Matriz_t *U, real_t *y,double *tx){

  int n=U->n; //tamanho da matriz U - > decomposição lU

  real_t *x = (real_t *) malloc(n * sizeof(real_t)); //vetor de solução X
  real_t soma;

  n--;
  x[n] = y[n]/U->A[n][n]; // y retornado da operação LY=B
  for(int linha = n-1; linha >=0; linha--){
    soma = 0;
    for(int coluna = linha; coluna <= n ;coluna++){
      soma += U->A[linha][coluna]*x[coluna];
    }
    x[linha] = (y[linha]-(soma))/U->A[linha][linha]; //ver se pode ter um zero
  }

  return x;
 
}//OK


//####################CÁLCULO DA INVERSA DADO L E U DA DECOMPOSIÇÃO LU####################
 int inversa(Matriz_t *L, Matriz_t *U,int n,double *tTotal, Matriz_t *ident){

  int linha, coluna;
  double time, elapsed;
  real_t  soma;
  real_t *y = (real_t *) malloc(n * sizeof(real_t));
  real_t *x = (real_t *) malloc(n * sizeof(real_t));

  
//-------------------------------------------------------------------------
  Matriz_t *inversa = alocaMatriz(n); //aloca a matriz que vai receber a inversa
  int j = 0;

  for (int col=0; col<n;col++){   // for na matriz identidade, para o vetor x receber a coluna n para o cálculo de (LYn=B)
    for (int lin = 0; lin<n; lin++){
      x[lin] = ident->A[lin][col];
    }


    time = timestamp();
    y = triangInferior(L,x,tTotal); //chama a função para a solução de (LY=B)
    elapsed = timestamp() - time;
    tTotal[0] += elapsed;


    time = timestamp();
    x = triangSuperior(U,y,tTotal); // chama a função para a solução de (UX=Y) passando y, calculado acima
    elapsed = timestamp() - time;
    tTotal[1] += elapsed;


    for (int i=0;i<n;i++){
      inversa->A[i][j]=x[i];  // a inversa recebe a coluna n da solução (UXn=Y)
    }
    j++;
  }

	for (int i=0;i<n;i++){
		for(int j=0;j<n;j++){
    		U->A [i][j] = inversa->A[i][j];  //retorna a inversa na MAIN
		}
  	}

  	free(ident);
  	
}//ok


//####################DECOMPOSIÇÃO DE UMA MATRIZ A EM LU (A=LU)####################
int decomposicaoLU(Matriz_t *matriz, double *tTotal){
	
	unsigned int n = matriz->n;
	int coluna, linha;
	Matriz_t *L = alocaMatriz(n);
	L->n = n;

	
	Matriz_t *ident = alocaMatriz(n);

  	for (linha=0; linha<n; linha++){
    	for (coluna = 0; coluna<n; coluna++){
      		if(linha == coluna){
        		ident->A[linha][coluna] = 1;
      		}else{
        		ident->A[linha][coluna] = 0;
      		}
    	}
  	}

	for (coluna = 0; coluna < n; coluna++){
		L->A[coluna][coluna] = 1; //fusão
		for (linha = coluna+1; linha < n; linha++){
			L->A[linha][coluna] = matriz->A[linha][coluna] / matriz->A[coluna][coluna]; //matriz L que recebe os coeficientes, Axx/pivô
			matriz->A[linha][coluna] = 0.0; //zera o coef na matriz U

			for (int j=coluna+1;j<n;j++){
				matriz->A[linha][j] = matriz->A[linha][j]-L->A[linha][coluna]*matriz->A[coluna][j]; //faz as operações no restante das colunas
			}
		}
	}

//fusão no for anterior
//  for (int lin = 0; lin < matriz->n; lin++){
//  	for (int col=0; col < matriz->n ; col++){
//      if(lin==col)
//        L->A[lin][col]=1;  //coloca 1 na diagonal principal de L
//    }
//  }

  //determinante recebe 1 pq não temos pivoteamento
  int determinante =1;

  determinante = calculaDeterminante(matriz,determinante);
  if(determinante==0){
  	printf("A matriz não possui inversa\n");
  }else
  	printf("\n");
    inversa(L,matriz,n,tTotal,ident);  //chama a função para cálculo da inversa

}//ok


//####################TROCA AS LINHA NO PIVOTEAMENTO PARCIAL####################
//NÃO UTILIZA MAIS
Matriz_t* trocaLinhas(Matriz_t *matriz, int linha_pivo, int coluna){

  double aux, tmp;

  for(int i =0; i< matriz->n; i++){
    aux = matriz->A[linha_pivo][i];
    matriz->A[linha_pivo][i] = matriz->A[coluna][i];
    matriz->A[coluna][i] = aux;
  }

  return matriz;
}//OK


//####################ENCONTRA O PIVÔ DA COLUNA PARA O PIVOTEAMENTO PARCIAL####################
int encontraMax(Matriz_t *matriz, int n,int linha, int coluna){

  real_t max;
  max = matriz->A[linha][coluna];
  int linha_pivo = 0 ;

  for(int i = linha; i < n-1; ++i) {
    if (fabs(matriz->A[i][coluna]) >= fabs(max)){
      max = matriz->A[i][coluna];
      linha_pivo = i;
    }
  }
  return linha_pivo;
}//OK


//####DECOMPOSIÇÃO DE UMA MATRIZ A EM LU (A=LU) UTILIZANDO A TÉCNICA DE PIVOTEAMENTO PARCIAL########
//Feito o pivoteamento com ponteiro
int decomposicaoLUpivoteamento(Matriz_t *matriz, double *tTotal){

  unsigned int n = matriz->n;
  int linha_pivo=0, coluna, linha,i ;
  int determinante = 1;
  
  real_t *aux;
  Matriz_t *L = alocaMatriz(n);
  L->n = n;


  Matriz_t *ident = alocaMatriz(n);

  for (linha=0; linha<n; linha++){
    for (coluna = 0; coluna<n; coluna++){
      if(linha == coluna){
        ident->A[linha][coluna] = 1;
      }else{
        ident->A[linha][coluna] = 0;
      }
    }
  }


  //aquiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiiii
  for (coluna = 0; coluna < n; coluna++){
    linha_pivo = encontraMax(matriz,n,coluna,coluna); //procura o maior valor para ser o pivô da coluna
    if (linha_pivo != coluna){ //quando linha e coluna são iguais é pq é a diag principal então não precisa trocar
      
//------OTIMIZAÇÃO ao envés de chamar a troca_linha, faz a troca de ponteiro sem chamar outra função--------------------------------
    	aux = matriz->A[coluna];
    	matriz->A[coluna] = matriz->A[linha_pivo];
    	matriz->A[linha_pivo] = aux;
    	//trocaLinhas(matriz, linha_pivo, coluna);  //troca as linhas quando o maior valor não está na coluna de cálculo
    	determinante = -1; // se troca linha inverte o sinal do determinante

    	//pivotear a identidade
   		aux = ident->A[coluna];
    	ident->A[coluna] = ident->A[linha_pivo];
    	ident->A[linha_pivo] = aux;
    }

    L->A[coluna][coluna] = 1; //fusão
	for (linha = coluna+1; linha < n; linha++){
		L->A[linha][coluna] = matriz->A[linha][coluna] / matriz->A[coluna][coluna]; //matriz L que recebe os coeficientes, Axx/pivô
		matriz->A[linha][coluna] = 0.0; //zera o coef na matriz U

		for (int j=coluna+1;j<n;j++){
			matriz->A[linha][j] = matriz->A[linha][j]-L->A[linha][coluna]*matriz->A[coluna][j]; //faz as operações no restante das colunas
		}
	}
  }

  determinante = calculaDeterminante(matriz,determinante);
  if(determinante==0){
  	printf("A matriz não possui inversa\n");
  	exit(0);
  }else
  	inversa(L,matriz,n,tTotal,ident); //chama o cálculo da inversa passando L U(matriz) e o tamanho
}//precisa arrumar