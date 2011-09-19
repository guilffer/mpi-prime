#include <stdio.h>
#include <string.h>
#include <mpi.h>
#include <gmp.h>
#include <stdlib.h>
#define TRUE 1
#define FALSE 0
#define INTERVALO 1000 //Valor do intervalo utilizado pelos nos para procurar por numeros primos
#define forever for(;;)

int verificaPrimo(int inf);

int main(int argc, char *argv[]) {

  int ret, rank, size,tag,i, j;
  int pos1, resultado;
  int superior, inferior, anterior, anterior2;
  int aux;
  int menosUm = -1;
  int atual;
  int contador;

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
  MPI_Status status;
  MPI_Request request;

  FILE * arq;

  ret = MPI_Init(&argc, &argv);

  ret = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  ret = MPI_Comm_size(MPI_COMM_WORLD, &size);

  tag=100;

  //Arquivo onde sao salvos os numeros primos encontrados
  arq = fopen("encontrados.txt", "w+");

  //Se o rank for igual a zero indica que o nó é o mestre.
  if (rank  == 0) {
    double t1, t2; //variaveis de tempo

    /*
     * Este for é responsavel por enviar os limites inferiores dos 
     * intervalos para cada um dos nós
     */
    for (i=0; i<=7; i++) {
      j = i * INTERVALO + 1;
      ret = MPI_Send( &j, 1, MPI_INT, i+1, tag, MPI_COMM_WORLD);
      limites[i+1] = j + INTERVALO;
      printf("\nIntervalo de %d e %d", i, limites[i+1]);
    }

    j = j + INTERVALO; //Armazena o ultimo intervalo enviado 

    //Envia uma mensagem de "start" para cada um dos nos ativos
    for( i=0; i<= 7 ; i++)
      ret = MPI_Send(&menosUm, 1, MPI_INT, i+1, tag, MPI_COMM_WORLD);
    
    //Armazena o valor do tempo quando comecou a busca paralela pelos numeros primos
    t1 = MPI_Wtime(); 

    /*
     * Este loop infinito possibilita o recebimento das mensagens enviadas pelos 
     * nos quando um destes acha um numero primo
     */
    forever {

      //Recebe o numero primo encontrado e armazena na variavel pos1
      ret = MPI_Recv(&pos1, 1, MPI_INT, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);

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
       * menor que o numero encontrado o no mestre envia para o host o novo valor do intervalo
       * a ser atualizado.
       */
      for (i=1; i<8; i++) {
        if (limites[i] < pos1) {
          //Envia o novo intervalo para o host que possui limite menor que o numero primo encontrado
          ret = MPI_Send(&j, 1, MPI_INT, i, tag, MPI_COMM_WORLD);

          //Atualiza o intervalo
          j = j+INTERVALO;

          //Atualiza o limite do no que foi atualizado
          limites[i] = j;
        }
      }
    }
  
  //Caso nao seja o no mestre o no entrara no else abaixo
  } else {
    //Recebe o limite inferior inicial
    ret = MPI_Recv(&inferior, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);

    //Calcula o limite superior
    superior = inferior + INTERVALO;

    /*
     * Armazena o valor do limite inferior para posteriormente 
     * descobrir se o intervalo deve ser atualizado ou nao
     */
    anterior = anterior2 = inferior;

    //Recebe a mensagem de "start"
    ret = MPI_Recv(&aux, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &status);

    //Inicia a busca pelo numero primo
    forever {
      //Enquanto o limite inferior do intervalo for menor que o limite superior     
      while (inferior < superior) {
        //Chama a funcao que busca pelo primo 
        resultado = verificaPrimo(inferior);

        //Se o resultado diferente de -1 o numero é primo
        if (resultado != -1) {
          //Envia o numero primo para o no mestre
          ret = MPI_Send(&resultado, 1, MPI_INT, 0, tag, MPI_COMM_WORLD);
        }

        /*
         * Atualiza o limite inferior. 
         * Vale lembrar que pela regra do numero primo de Mersenne o expoente da
         * regra (2^n)-1 é sempre impar, sendo assim o limite inferior 
         * é atualizado de dois em dois 
         */
        inferior +=2;
        
        //Realiza o receive nao bloqueante, que verifica se o no mestre enviou um novo intervalo ou nao
        ret = MPI_Irecv(&anterior, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &request);

        /*
         * Se o valor recebido pelo receive nao bloqueante for diferente do limite inferior 
         * recebido anteriormente o intervalo deve ser atualizado
         */
        if (anterior != anterior2) {
          anterior2 = anterior;
          inferior = anterior; //Atualiza o limite inferior 
          superior = inferior + INTERVALO; //Atualiza o limite superior
          printf("\nNo %d - Atualizei o intervalo para %d", rank, inferior);
          fflush(stdout);
        }
      }

      /*
       * Caso o intervalo tenha terminado o nó fica esperando por um novo intervalo, a ser enviado
       * pelo no mestre
       */
      ret = MPI_Irecv(&anterior, 1, MPI_INT, 0, tag, MPI_COMM_WORLD, &request);
      if (anterior != anterior2) {
        anterior2 = anterior;
        inferior = anterior;
        superior = inferior + INTERVALO;
        printf("\nNo %d - Atualizei o intervalo para %d", rank, inferior);
        fflush(stdout);
      }
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
