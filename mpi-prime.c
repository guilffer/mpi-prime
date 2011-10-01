#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <gmp.h>
#include <stdlib.h>
#define TRUE 1
#define FALSE 0
#define INTERVALO 100000
#define forever for(;;)

#define mpi_send(X,NODE,MPI) MPI_Send(&X, 1, MPI_INT, NODE, MPI.tag, MPI_COMM_WORLD)
#define mpi_recv(X,NODE,MPI) MPI_Recv(&X, 1, MPI_INT, NODE, MPI.tag, MPI_COMM_WORLD, &MPI.status)
#define mpi_recv_any(X,MPI) MPI_Recv(&X, 1, MPI_INT, MPI_ANY_SOURCE, MPI.tag, MPI_COMM_WORLD, &MPI.status)
#define mpi_irecv(X,NODE,MPI) MPI_Irecv(&X, 1, MPI_INT, NODE, MPI.tag, MPI_COMM_WORLD, &MPI.request)

typedef struct {
  int superior, inferior, anterior;
} intervalo;

typedef struct {
  int ret, rank, size, tag;
  MPI_Status status;
  MPI_Request request;
} mpi_data;

int verificaPrimo(int inf);
void recebeNovoIntervalo(intervalo *it);

int main(int argc, char *argv[]) {

  int i, j;
  int pos1, resultado;
  int dump;
  int startSignal = -1;
  int contador;
  intervalo it;
  mpi_data mpi;
  FILE * arq;

  /*
   * Vetor de inteiros que armazena os limites de cada um dos hosts para que 
   * o no mestre possa realizar o balanceamento de carga automatico para cada
   * um dos nós.
   */
  long int limites[10];

  //variaveis da biblioteca gmp utilizadas durante a execucao
  mpz_t numeroGeral, numeroCorrente;
  //Inicializacao das variaveis
  mpz_init(numeroGeral);
  mpz_init(numeroCorrente);

  //Variaveis de status e request utilizadas pelo MPI

  mpi.ret = MPI_Init(&argc, &argv);
  mpi.ret = MPI_Comm_rank(MPI_COMM_WORLD, &mpi.rank);
  mpi.ret = MPI_Comm_size(MPI_COMM_WORLD, &mpi.size);

  arq = fopen("encontrados.txt", "w+");

  //Se o rank for igual a zero indica que o nó é o mestre.
  if (mpi.rank  == 0) {
    double t1, t2; //variaveis de tempo

    /*
     * Este for é responsavel por enviar os limites inferiores dos 
     * INTERVALOs para cada um dos nós
     */
    for (i=0; i<mpi.size; i++) {
      j = i * INTERVALO + 1;
      mpi.ret = mpi_send(j, i+1, mpi);
      limites[i+1] = j + INTERVALO;
      printf("\nIntervalo de %d e %d", i, limites[i+1]);
    }

    j = j + INTERVALO; //Armazena o ultimo INTERVALO enviado 

    //Envia uma mensagem de "start" para cada um dos nos ativos
    for(i=0; i<mpi.size; i++)
      mpi.ret = mpi_send(startSignal, i+1, mpi);
    
    //Armazena o valor do tempo quando comecou a busca paralela pelos numeros primos
    t1 = MPI_Wtime(); 

    /*
     * Este loop infinito possibilita o recebimento das mensagens enviadas pelos 
     * nos quando um destes acha um numero primo
     */
    forever {

      //Recebe o numero primo encontrado e armazena na variavel pos1
      mpi.ret = mpi_recv_any(pos1, mpi);

      //Pega o tempo do recebimento
      t2 = MPI_Wtime();

      //Imprime o numero primo de Mersenne e apresenta o tempo em que foi encontrado
      printf("\nO numero M%d e primo - Tempo = %f",pos1, t2-t1);
      fflush(stdout);

      //Escreve o numero encontrado e o valor do tempo no arquivo
      fprintf(arq, "\nO numero M%d e primo - Tempo = %f",pos1, t2-t1);
      fflush(arq);

      /*
       * Este for é responsavel por verificar os limites superiores de cada um dos nós 
       * e realizar o balanceamento de carga automatico. Caso um dos limites dos nos seja 
       * menor que o numero encontrado o no mestre envia para o host o novo valor do INTERVALO
       * a ser atualizado.
       */
      for (i=1; i<mpi.size; i++) {
        if (limites[i] < pos1) {
          //Envia o novo INTERVALO para o host que possui limite menor que o numero primo encontrado
          mpi.ret = mpi_send(j, i, mpi);
          //Atualiza o INTERVALO
          j = j+INTERVALO;
          //Atualiza o limite do no que foi atualizado
          limites[i] = j;
        }
      }
    }
  
  //Caso nao seja o no mestre o no entrara no else abaixo
  } else {
    //Recebe o limite inferior inicial
    mpi.ret = mpi_recv(it.inferior, 0, mpi);

    //Calcula o limite superior
    it.superior = it.inferior + INTERVALO;

    /*
     * Armazena o valor do limite inferior para posteriormente 
     * descobrir se o INTERVALO deve ser atualizado ou nao
     */
    it.anterior = it.inferior;

    //Recebe a mensagem de "start"
    mpi.ret = mpi_recv(dump, 0, mpi);

    //Inicia a busca pelo numero primo
    forever {
      //Enquanto o limite inferior do INTERVALO for menor que o limite superior     
      while (it.inferior < it.superior) {
        //Chama a funcao que busca pelo primo 
        resultado = verificaPrimo(it.inferior);

        //Se o resultado diferente de -1 o numero é primo
        if (resultado != -1) {
          //Envia o numero primo para o no mestre
          mpi.ret = mpi_send(resultado, 0, mpi);
        }

        /*
         * Atualiza o limite inferior. 
         * Vale lembrar que pela regra do numero primo de Mersenne o expoente da
         * regra (2^n)-1 é sempre impar, sendo assim o limite inferior 
         * é atualizado de dois em dois 
         */
        it.inferior +=2;
        
        recebeNovoIntervalo(&it);
      }

      /*
       * Caso o INTERVALO tenha terminado o nó fica esperando por um novo INTERVALO, a ser enviado
       * pelo no mestre
       */
      recebeNovoIntervalo(&it);
    }
  }

  fclose(arq);
  return 0;
}


/*
 * Funcao que realiza a busca pelo numero primo de Mersenne 
 * utilizando o método de Lucas Lehmer
 */
int verificaPrimo(int inf){
  int p, i;
  mpz_t lol, aux, pop;

  mpz_init(lol);
  mpz_init(aux);
  mpz_init(pop);

  p = inf;

  mpz_set_ui(lol,(long)(4));
  mpz_ui_pow_ui(aux,(long)(2),p);
  mpz_sub_ui(pop,aux,(long)(1));

  for (i=1; i+2<=p; i++) {
    mpz_powm_ui(aux,lol,(long)(2),pop);
    mpz_sub_ui(lol,aux,(long)(2));
  }

  return (mpz_sgn(lol) == 0) ? p : -1;
}

void recebeNovoIntervalo(intervalo *it) {
  int aux = it->anterior;
  mpi.ret = mpi_irecv(it->anterior, 0, mpi);
  if (it->anterior != aux) {
    aux = it->anterior;
    it->inferior = it->anterior; //Atualiza o limite inferior 
    it->superior = it->inferior + INTERVALO; //Atualiza o limite superior
    printf("\nNo %d - Atualizei o INTERVALO para %d", mpi.rank, it->inferior);
    fflush(stdout);
  }
}
