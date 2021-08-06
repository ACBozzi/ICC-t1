//Além disso, ao final de cada bloco (após a matriz inversa) deve constar:
//O tempo de execução da triangularização;
//O tempo de cálculo de Y e X nas fases Ly=b e Ux=y
//O valor da norma L2 do resíduo de A.A' = I, isto é, a norma L2 do resíduo de cada coluna da matriz 
//identidade com a coluna correspondente da matriz A.A', onde A'  é o resultado do cálculo de A-1 pelo Método da Fatoração LU.

//Desenvolvolvedores:
	//Anna Caroline Bozzi - acb17 -  GRR20173532
	//João Vitor Gualberto Figueiredo - jvgf17 - GRR20177037

#ifndef __SISLINEAR_H__
#define __SISLINEAR_H__

typedef float real_t;

typedef struct {
  unsigned int n; // tamanho da matriz
  real_t **A; // coeficientes
} Matriz_t;

//Alocaçao  de memória
Matriz_t* alocaMatriz (unsigned int n);
void copiaMatriz(Matriz_t *SL, Matriz_t *original);

//Leitura e impressão de sistemas lineares
Matriz_t *lerMatriz ();
void prnMatriz (Matriz_t *matriz, int file_write, FILE *file);
void prnVetor (real_t *vet, unsigned int n);

int calculaDeterminante(Matriz_t *matriz, int determ);
real_t* normaL2residuo(Matriz_t *original, Matriz_t *inversa, real_t *resi);

//cálculo da decomṕsição de A em LU
int decomposicaoLU (Matriz_t *matriz, double *tTotal);

//cálculo da inversa dados L e U
int inversa(Matriz_t *L, Matriz_t *matriz, int n,double *tTotal, Matriz_t *identidade);

//solução de ly=b
real_t* triangSuperior(Matriz_t *matriz, real_t *y,double *tx);
//solução de ux=y
real_t* triangInferior(Matriz_t *L, real_t *x,double *ty);

//cálculo da fatoração com pivoteamento parcial, busca de máximo e troca de linhass
int encontraMax(Matriz_t *matriz, int n,int linha, int coluna);
Matriz_t* trocaLinhas(Matriz_t *matriz, int linha_pivo, int coluna);
int decomposicaoLUpivoteamento(Matriz_t *matriz, double *tTotal);

//real_t normaL2residuo(Matriz_t *original, Matriz_t *inversa);

#endif // __SISLINEAR_H__