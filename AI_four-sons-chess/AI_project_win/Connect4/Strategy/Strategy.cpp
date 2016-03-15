#include <iostream>
#include <memory>
#include "Point.h"
#include "Strategy.h"

/*#define dead1 0;
#define alive1 2;
#define dead2 3;
#define alive2 100;
#define dead3 500000;
#define alive3 500000;//��ֵĸ���״̬��Ȩֵ����1��ֻ��һ���������������Ա߻���һ���Է����Ӷ��ţ���2����������������������û�жԷ����Ӷ��š���
*/
const int dead1 = 10;
const int alive1 = 30;
const int dead2 = 50;
const int alive2 = 8000;
const int dead3 = 9000;
const int alive3 = 100000;
const int dead4 = 20000000;//���һ��
const int alive4 = 20000000;
const int depth = 5;//alpha_beta�㷨��Ѱ���ĸ߶ȣ����
//const int legth = 50000;//(int)((pow(12.0, 5*1.0 + 1) - 1) / 11);//�ڵ������Ŀ
using namespace std;

int Y = -1;
bool alphaB = true;
struct node
{
public:
    node *father;
	//node *sons[12];//�����12������
	int alpha;//�������½�
	int beta;//��С����Ͻ�
	bool isMax;//�����С
   // bool hasfp;//Ȩֵ
	int depth;//���
	int **b;//�õ������״̬
	//int x, y;//����õ����ı����������
public:
	node()
	{
		father = NULL;
		//for(int i = 0; i < 12; i++)
		//{
			//sons[i] = NULL;
		//}
		alpha = -100000000;
		beta = 1000000000;
		isMax = true;
		//hasfp = false;//Ҷ�ڵ����fp�������ڵ�ֻ��������½磬Ҳ�������ж��Ƿ��Ѽ�������½磬Ҳ�����ж������Ǵ������ݻ�������
		depth = 0;
		//x = 0;
		//y = 0;
		b = new int*[12];
		for(int i = 0; i < 12; i++)
		{
			b[i] = new int[12];
			for(int j = 0; j < 12; j++)
				b[i][j] = 0;
		}
	}
	~node()
	{
		for(int i = 0; i < 12; i++)
		{
			delete[]b[i];

		}
		delete []b;
	}
};



int *getND(int **b, int max, int M, int N, int noX, int noY)//M N ����ֳ��ȺͿ��
{
	//cout<<"M  :  "<<M<<"   "<<"N  :  "<<N<<endl;
	int i, j;
	int *num;
	num = new int[8];
	memset(num, 0, sizeof(num));
	/*for(int m = 0; m < M; m++)
	{
		cout<<"for   "<<endl;
		for(int n = 0; n < N; n++)
		{
			cout<<"b...."<<" ";
			cout<<b[m][n]<<" ";
		}
		cout<<endl;
	}
	cout<<"after cout board  **************"<<endl;*/
	if(max == 1)//�����ҷ���״̬��
	{
		//cout<<"max == 1  "<<endl;
		for(i = 0; i < M; i++)//���ţ�һ��һ�м���
		{
			for(j = 0; j < N; j++)
			{
				if(b[i][j] == 2)//�����ҷ�����
				{
					if(j == N - 1)
						num[0]++;//��һ
					else
					{
						if(j == 0)//��
						{
							int k;
							for(k = j + 1; k < N; k++)
							{
								if(b[i][k] != 2)
									break;
							}
							num[2*(k - j) - 2]++;
						}
						else
						{
							int k;
							for(k = j + 1; k < N; k++)
							{
								if(b[i][k] != 2)
									break;
							}
							if(b[i][j - 1] == 1 || (i == noX && j - 1 == noY))//�Է����ӻ�ô��������Ӷ��ţ�Ϊ��,,,,�˴������˲������ӵ�
							{
								num[2*(k - j) - 2]++;
							}
							else
							{
								if((k < N) && (b[i][k] == 0) && !(i == noX && k == noY))//�ұ߽�Ҳ���ţ���
									num[2*(k - j) - 1]++;
								else
									num[2*(k - j) - 2]++;//��
							}
						}
					}
				}
			}
		}
		//cout<<"after һ��һ�У�����"<<endl;
		/*for(int t = 0; t < 10; t++)
		{
			cout<<num[t]<<endl;
		}*/
		for(j = 0; j < N; j++)
		{
			for(i = 0; i < M; i++)//һ��һ��
			{
				if(b[i][j] == 2)//�����ҷ�����
				{
					if(i == M - 1)
						num[0]++;//��һ
					else
					{
						if(i == 0)//��
						{
							int k;
							for(k = i + 1; k < M; k++)
							{
								if(b[k][j] != 2)
									break;
							}
							num[2*(k - i) - 2]++;
						}
						else
						{
							int k;
							for(k = i + 1; k < M; k++)
							{
								if(b[k][j] != 2)
									break;
							}
							if(b[i - 1][j] == 1 || (i - 1 == noX && j == noY))//�Է����ӻ�ô��������Ӷ��ţ�Ϊ��,,,,�˴������˲������ӵ�
							{
								num[2*(k - i) - 2]++;
							}
							else
							{
								if((k < M) && (b[k][j] == 0) && !(k == noX && j == noY))//�ұ߽�Ҳ���ţ���
									num[2*(k - i) - 1]++;
								else
									num[2*(k - i) - 2]++;//��
							}
						}
					}
				}
			}
		}
		for(i = 0; i < M; i++)//б�ţ� ��һ���ҵ�����һ����ĵ�һ������һ��һ������
		{
			for(j = 0; j < N; j++)
			{
				if(b[i][j] == 2)
				{
					int k;
					if(i == M - 1)
						break;
					else
					{
						for(k = 0; k + j < N && k + i < M; k++)//�����Ͻǵ����½�
						{
							if(b[i + k][j + k] != 2)
							{
								break;
							}
							if((i > 0 && j > 0) && b[i - 1][j - 1] == 0 && !(i - 1 == noX && j - 1 == noY) && i + k < M && j + k < N && b[i + k][j + k] == 0 && !(i + k == noX && j + k == noY))//�������
							{
								//if(k > 1)
									num[2 * k - 1]++;
							}
							else//��
							{
								//if(k > 1)//ȥ������һ�����ӵ����
									num[2*k - 2]++;
							}
						}
						for(k = 0; i + k < M && j - k >= 0; k++)
						{
							if(b[i + k][j - k] != 2)
								break;
							if(i > 0 && j < N - 1 && b[i - 1][j + 1] == 0 && !(i - 1 == noX && j + 1 == noY) && i + k < M && j - k >= 0 && b[i + k][j - k] == 0 && !(i + k == noX && j - k == noY))//�������
							{
								//if(k > 1)
									num[2*k - 1]++;
							}
							else
							{
								//if(k > 1)
									num[2*k - 2]++;
							}
						}
					}
				}
			}
		}
	}
	else
	{
		for(i = 0; i < M; i++)//���ţ�һ��һ�м���
		{
			for(j = 0; j < N; j++)
			{
				if(b[i][j] == 1)//�����ҷ�����
				{
					if(j == N - 1)
						num[0]++;//��һ
					else
					{
						if(j == 0)//��
						{
							int k;
							for(k = j + 1; k < N; k++)
							{
								if(b[i][k] != 1)
									break;
							}
							num[2*(k - j) - 2]++;
						}
						else
						{
							int k;
							for(k = j + 1; k < N; k++)
							{
								if(b[i][k] != 1)
									break;
							}
							if(b[i][j - 1] == 2 || (i == noX && j - 1 == noY))//�Է����ӻ�ô��������Ӷ��ţ�Ϊ��,,,,�˴������˲������ӵ�
							{
								num[2*(k - j) - 2]++;
							}
							else
							{
								if((k < N) && (b[i][k] == 0) && !(i == noX && k == noY))//�ұ߽�Ҳ���ţ���
									num[2*(k - j) - 1]++;
								else
									num[2*(k - j) - 2]++;//��
							}
						}
					}
				}
			}
		}
		for(j = 0; j < N; j++)
		{
			for(i = 0; i < M; i++)//һ��һ��
			{
				if(b[i][j] == 1)//�����ҷ�����
				{
					if(i == M - 1)
						num[0]++;//��һ
					else
					{
						if(i == 0)//��
						{
							int k;
							for(k = i + 1; k < M; k++)
							{
								if(b[k][j] != 1)
									break;
							}
							num[2*(k - i) - 2]++;
						}
						else
						{
							int k;
							for(k = i + 1; k < M; k++)
							{
								if(b[k][j] != 1)
									break;
							}
							if(b[i - 1][j] == 2 || (i - 1 == noX && j == noY))//�Է����ӻ�ô��������Ӷ��ţ�Ϊ��,,,,�˴������˲������ӵ�
							{
								num[2*(k - i) - 2]++;
							}
							else
							{
								if((k < M) && (b[k][j] == 0) && !(k == noX && j == noY))//�ұ߽�Ҳ���ţ���
									num[2*(k - i) - 1]++;
								else
									num[2*(k - i) - 2]++;//��
							}
						}
					}
				}
			}
		}
		for(i = 0; i < M; i++)//б�ţ� ��һ���ҵ�����һ����ĵ�һ������һ��һ������
		{
			for(j = 0; j < N; j++)
			{
				if(b[i][j] == 1)
				{
					int k;
					if(i == M - 1)
						break;
					else
					{
						for(k = 0; k + j < N && k + i < M; k++)//�����Ͻǵ����½�
						{
							if(b[i + k][j + k] != 1)
							{
								break;
							}
							if((i > 0 && j > 0) && b[i - 1][j - 1] == 0 && !(i - 1 == noX && j - 1 == noY) && i + k < M && j + k < N && b[i + k][j + k] == 0 && !(i + k == noX && j + k == noY))//�������
							{
								//if(k > 1)
									num[2 * k - 1]++;
							}
							else//��
							{
								//if(k > 1)//ȥ������һ�����ӵ����
									num[2*k - 2]++;
							}
						}
						for(k = 0; i + k < M && j - k >= 0; k++)
						{
							if(b[i + k][j - k] != 1)
								break;
							if(i > 0 && j < N - 1 && b[i - 1][j + 1] == 0 && !(i - 1 == noX && j + 1 == noY) && i + k < M && j - k >= 0 && b[i + k][j - k] == 0 && !(i + k == noX && j - k == noY))//�������
							{
								//if(k > 1)
									num[2*k - 1]++;
							}
							else
							{
								//if(k > 1)
									num[2*k - 2]++;
							}
						}
					}
				}
			}
		}
	}
	return num;
}

int getfp(int **b, int M, int N, int noX, int noY)//�ú������ݵ�ǰ���״̬�����״̬�µ�����ֵfp = max - min
{
	/*for(int m = 0; m < M; m++)
	{
		for(int n = 0; n < N; n++)
		{
			cout<<b[m][n]<<" ";
		}
		cout<<endl;
	}*///OK!!!!
	int *numOfADmax;//����״̬��������״̬��1�±���0��״̬��1�±���1��������
	numOfADmax = new int[8];
	memset(numOfADmax, 0, sizeof(numOfADmax));
	//cout<<"getND   "<<endl;
	numOfADmax = getND(b, 1, M, N, noX, noY);//��/����/max ����״̬����
	//cout<<"�ҷ���״̬������֪"<<endl;
	//for(int t = 0; t < 8; t++)
		//cout<<numOfADmax[t]<<" ";
	//cout<<"after cout numOfADmax "<<endl;
	int *numOfADmin;
	numOfADmin = new int[8];
	numOfADmin = getND(b, 0, M, N, noX, noY);//�Է���״̬����
   // for(int t = 0; t < 8; t++)
	//	cout<<numOfADmin[t]<<" ";
	//cout<<"after cout numOfADmin "<<endl;
	int fp, i, j;
	int max = 0;
	int min = 0;
	//cout<<"alive1 "<<alive1 << numOfADmin[1] * alive1<< endl;
	max = max + numOfADmax[0] * dead1 + numOfADmax[1] * alive1 + numOfADmax[2] * dead2 + numOfADmax[3] * alive2 + numOfADmax[4] * dead3 + numOfADmax[5] * alive3 + numOfADmax[6] * dead4 + numOfADmax[7] * alive4;
	min = min + numOfADmin[0] * dead1 + numOfADmin[1] * alive1 + numOfADmin[2] * dead2 + numOfADmin[3] * alive2 + numOfADmin[4] * dead3 + numOfADmin[5] * alive3 + numOfADmin[6] * dead4 + numOfADmin[7] * alive4;
	//if(min > 10000000)
		//max = 0;
	fp = max - min;
	//cout<<max <<" "<<min<<endl;
	return fp;
}

int alphaBeta(node *root, int M, int N, int noX, int noY, int *_top)
{
	//cout<<"in to alphaBeta   .............."<<endl;
	if(root->depth == 0)//Ҷ�ڵ�
		return getfp(root->b, M, N, noX, noY);
	int i, j;
	int x, y;
	int value;//��Ҫ�����ص�ֵ��
	//cout<<"root :    "<<endl;
	//cout<<root->depth<<"  "<<root->isMax<<endl;
	for(j = 0; j < N; j++)
	{
		if((_top[j] - 1 >= 0) && (_top[j] - 1 < M) && !(_top[j] - 1 == noX && j == noY))//���п�������
		{
			x = _top[j] - 1;
			y = j;
			//cout<<"x,    y   "<<x<<"  "<<y<<endl;
			_top[j]--;
			node _next;
			node *next = &_next;
			for(int m = 0; m < M; m++)
			{
				for(int n = 0; n < N; n++)
					next->b[m][n] = root->b[m][n];
			}
			//root->sons[j] = next;
			next->father = root;
			//next->x = x;
			//next->y = y;
			next->depth = root->depth - 1;//����ڼ���
			//cout<<"root ....."<<endl;
			if(root->isMax)//���ڵ�Ϊ����ڵ�
			{
				//cout<<"isMax  ...."<<endl;
				next->isMax = false;

				next->b[x][y] = 2;//����ڵ����ӱ�ɼ�С�ڵ�
				if(j == 0)//����ߵĺ���
				{
					/*if(getfp(next->b, M, N, noX, noY) > 1000000000 && root->depth == depth)//��ʤ
					{	
						Y = j;
						cout<<"getfp(next->b, M, N, noX, noY)   "<<getfp(next->b, M, N, noX, noY)<<endl;
						break;
					}*/
					//cout<<"�����.....  "<<endl;
					next->beta = alphaBeta(next, M, N, noX, noY, _top);
					//cout<<"�ݹ�֮��next->beta    "<<next->beta<<endl;
					root->alpha = next->beta;
					if(root->depth == depth)
						Y = j;
					//cout<<"Y    ::   "<<Y<<endl;
					value = next->beta;
				}
				else
				{
					//cout<<"next is not the first child  ......."<<endl;
					/*if(getfp(next->b, M, N, noX, noY) > 1000000000 && root->depth == depth)//��ʤ
					{	
						Y = j;
						cout<<"��ʤ  ..."<<getfp(next->b, M, N, noX, noY)<<endl;
						break;
					}*/
					//node _temp;
					node *temp = root->father;
					//temp = root->father;
					while(temp)
					{
						if(root->alpha >= temp->beta)//beta��֦
						{
							break;//����ֹͣ
						}
						if(temp->father)
						   temp = temp->father->father;
						else  temp=NULL;
					}
					if(temp)//��֦��
					{
						//cout<<"��֦   ����������"<<endl;
						value = root->alpha;
						//continue;
						_top[j]++;
						break;
					}
					else//????û�з���
					{
						//cout<<"û�м�֦"<<endl;
						next->beta = alphaBeta(next, M, N, noX, noY, _top);//������չ
						//cout<<"next -> beta  .....  "<<next->beta<<endl;
						//root->father->beta = root->alpha;//���ϸ��� 
						if(root->depth == depth)
						{
							if(next->beta > root->alpha)
							   Y = j;
						}
						if(next->beta > root->alpha)
						{
							value = next->beta;
							root->alpha = value;
						}
						else
							value = root->alpha;
					}
				}

				
			}
			else
			{
				//cout<<"root is not max  ....."<<endl;
				//cout<<"  depth  "<<next->depth<<endl;
				next -> isMax = true;
				next->b[x][y] = 1;
				/*if(getfp(next->b, M, N, noX, noY) < -1000000 && root->depth == depth - 1)//����
					{	
						value = getfp(next->b, M, N, noX, noY);
						cout<<"���ˡ�������"<<endl;
						break;
					}*/
				if(j == 0)//����ߵĺ���
				{
					//cout<<"the first child  ......"<<endl;
					next->alpha = alphaBeta(next, M, N, noX, noY, _top);
					root->beta = next->alpha;
					if(next->depth == depth)
						Y = j;
					//cout<<"Y   ....."<<endl;
					value = next->alpha;
					//cout<<"next->alpha  .....  "<<next->alpha<<endl;
					
				}
				else
				{
					//cout<<"is not the first .. "<<endl;
					node *temp = root->father;
					//temp = root->father;
					while(temp)
					{
						if(root->beta <= temp->alpha)//beta��֦
						{
							break;//����ֹͣ
						}
						if(temp->father)
						  temp = temp->father->father;
						else
							temp = NULL;
					}
					if(temp)//��֦��
					{
						//cout<<"��֦  ����������  "<<endl;
						value = root->beta;
						//continue;
						_top[j]++;
						break;

					}
					else//??????
					{
						next->alpha = alphaBeta(next, M, N, noX, noY, _top);//������չ
						//cout<<"next alpha  ....  "<<next->alpha<<endl;
						/*if(root->depth == depth)
						{
							//root->father->alpha = root->beta;//���ϸ���
							if(next->beta > root->alpha)
							   Y = j;
						}*///impossible
						if(next->alpha < root->beta)
						{
								value = next->alpha;
								root->beta = value;
						}
						else
							value = root->beta;
						   

					}
				}

				
			}
		//cout<<"this node ..........  "<<endl;
		//cout<<"next->depth : "<<next->depth<<endl;
		//cout<<"next->isMax : "<<next->isMax<<endl;
		//cout<<"next->alpha : "<<next->alpha<<endl;
		//cout<<"next->beta :  "<<next->beta<<endl;
		//cout<<"x,   y:       "<<next->x<<" "<<next->y<<endl;
		_top[j]++;
		//cout<<"_top[j]++ ..... "<<endl;
		//cout<<endl;
		//��ԭ
		}
	}
	return value;
}
/*
	���Ժ����ӿ�,�ú������Կ�ƽ̨����,ÿ�δ��뵱ǰ״̬,Ҫ�����������ӵ�,�����ӵ������һ��������Ϸ��������ӵ�,��Ȼ�Կ�ƽ̨��ֱ����Ϊ��ĳ�������
	
	input:
		Ϊ�˷�ֹ�ԶԿ�ƽ̨ά����������ɸ��ģ����д���Ĳ�����Ϊconst����
		M, N : ���̴�С M - ���� N - ���� ����0��ʼ�ƣ� ���Ͻ�Ϊ����ԭ�㣬����x��ǣ�����y���
		top : ��ǰ����ÿһ���ж���ʵ��λ��. e.g. ��i��Ϊ��,��_top[i] == M, ��i������,��_top[i] == 0//����´��ڴ������ӣ�������_top[i] - 1��
		_board : ���̵�һά�����ʾ, Ϊ�˷���ʹ�ã��ڸú����տ�ʼ���������Ѿ�����ת��Ϊ�˶�ά����board
				��ֻ��ֱ��ʹ��board���ɣ����Ͻ�Ϊ����ԭ�㣬�����[0][0]��ʼ��(����[1][1])
				board[x][y]��ʾ��x�С���y�еĵ�(��0��ʼ��)
				board[x][y] == 0/1/2 �ֱ��Ӧ(x,y)�� ������/���û�����/�г������,�������ӵ㴦��ֵҲΪ0
		lastX, lastY : �Է���һ�����ӵ�λ��, ����ܲ���Ҫ�ò�����Ҳ������Ҫ�Ĳ������ǶԷ�һ����
				����λ�ã���ʱ��������Լ��ĳ����м�¼�Է������ಽ������λ�ã�����ȫȡ�������Լ��Ĳ���
		noX, noY : �����ϵĲ������ӵ�(ע:��ʵ���������top�Ѿ����㴦���˲������ӵ㣬Ҳ����˵���ĳһ��
				������ӵ�����ǡ�ǲ������ӵ㣬��ôUI�����еĴ�����Ѿ������е�topֵ�ֽ�����һ�μ�һ������
				��������Ĵ�����Ҳ���Ը�����ʹ��noX��noY��������������ȫ��Ϊtop������ǵ�ǰÿ�еĶ�������,
				��Ȼ�������ʹ��lastX,lastY�������п��ܾ�Ҫͬʱ����noX��noY��)
		���ϲ���ʵ���ϰ����˵�ǰ״̬(M N _top _board)�Լ���ʷ��Ϣ(lastX lastY),��Ҫ���ľ�������Щ��Ϣ�¸������������ǵ����ӵ�
	output:
		������ӵ�Point
*/


extern "C" __declspec(dllexport) Point* getPoint(const int M, const int N, const int* top, const int* _board, 
	const int lastX, const int lastY, const int noX, const int noY){
	/*
		��Ҫ������δ���
	*/
	int x = -1, y = -1;//���ս�������ӵ�浽x,y��
	int** board = new int*[M];
	for(int i = 0; i < M; i++){
		board[i] = new int[N];
		for(int j = 0; j < N; j++){
			board[i][j] = _board[i * N + j];
		}
	}
	
	/*
		�������Լ��Ĳ������������ӵ�,Ҳ���Ǹ�����Ĳ�����ɶ�x,y�ĸ�ֵ
		�ò��ֶԲ���ʹ��û�����ƣ�Ϊ�˷���ʵ�֣�����Զ����Լ��µ��ࡢ.h�ļ���.cpp�ļ�
	*/
	//Add your own code below
	/*
     //a naive example
	for (int i = N-1; i >= 0; i--) {
		if (top[i] > 0) {
			x = top[i] - 1;
			y = i;
			break;
		}//ֱ�Ӵ��ұ�һ������
	}*/
    //����С*****
	
	int i, j;
	int tempf = -1000000000;//��С�������
	int xx, yy;//��¼���tempfʱ������x, y
	for(i = 0; i < N; i++)//��
	{
		if(top[i] > 0)
		{
			x = top[i] - 1;
			y = i;
			board[x][y] = 2;
			int temp = 2000000000;//��ʼ���������С
			int x1 = -1;
			int y1 = -1;
			int tempx = -1;
			int tempy = -1;
			for(j = 0; j < N; j++)
			{
				if(j != i && top[j] > 0)
				{
				   x1 = top[j] - 1;
				   y1 = j;
				}
				if(j == i && top[j] > 1)
				{
					x1 = top[j] - 2;
					y1 = j;
				}
				if(y1 > -1)//�ڸ����ҵ������ӵ�
				{
					board[x1][y1] = 1;
					int fp;
					//int max = 0;
					//int min = 0;

					/*for(int m = 0; m < M; m++)
					{
						for(int n = 0; n < N; n++)
						{
							if(board[m][n] == 2)
							{
								if((m > 2 && board[m - 1][n] != 1 && board[m - 2][n] != 1 && board[m - 3][n] != 1) || (m < M - 1 && m > 1 && board[m + 1][n] != 1 && board[m - 1][n] != 1 && board[m - 2][n] != 1) || (m > 0 && m < M - 2 && board[m + 2][n] != 1 && board[m + 1][n] != 1 && board[m - 1][n] != 1) || (m < M - 3 && board[m + 3][n] != 1 && board[m + 2][n] != 1 && board[m + 1][n] != 1))//���µ���
									max++;
								if((m > 2 && board[m - 1][n] != 2 && board[m - 2][n] != 2 && board[m - 3][n] != 2) || (m > 1 && m < M - 1 && board[m + 1][n] != 2 && board[m - 1][n] != 2 && board[m - 2][n] != 2) || (m > 0 && m < M - 2 && board[m + 2][n] != 2 && board[m + 1][n] != 2 && board[m - 1][n] != 2) || (m < M - 3 && board[m + 3][n] != 2 && board[m + 2][n] != 2 && board[m + 1][n] != 2))//���µ���
									min--;
								if((board[m][n - 1] != 1 && board[m][n - 2] != 1 && board[m][n - 3] != 1 && n > 2) || (n > 1 && n < N - 1 && board[m][n + 1] != 1 && board[m][n - 1] != 1 && board[m][n - 2] != 1) || (n > 0 && n < N - 2 && board[m][n + 2] != 1 && board[m][n + 1] != 1 && board[m][n - 1] != 1) || (n < N - 3 && board[m][n + 3] != 1 && board[m][n + 2] != 1 && board[m][n + 1] != 1))//������
									max++;
								if((board[m][n - 1] != 2 && board[m][n - 2] != 2 && board[m][n - 3] != 2 && n > 2) || (n > 1 && n < N - 1 && board[m][n + 1] != 2 && board[m][n - 1] != 2 && board[m][n - 2] != 2) || (n > 0 && n < N - 2 && board[m][n + 2] != 2 && board[m][n + 1] != 2 && board[m][n - 1] != 2) || (n < N - 3 && board[m][n + 3] != 2 && board[m][n + 2] != 2 && board[m][n + 1] != 2))//�����ң�����
									min--;
								if((m > 2 && board[m - 1][n] == 2 && board[m - 2][n] == 2 && board[m - 3][n] == 2) || (m < M - 1 && m > 1 && board[m + 1][n] == 2 && board[m - 1][n] == 2 && board[m - 2][n] == 2) || (m > 0 && m < M - 2 && board[m + 2][n] == 2 && board[m + 1][n] == 2 && board[m - 1][n] == 2) || (m < M - 3 && board[m + 3][n] == 2 && board[m + 2][n] == 2 && board[m + 1][n] == 2))
									max = max + 10;
								if((m > 2 && board[m - 1][n] == 1 && board[m - 2][n] == 1 && board[m - 3][n] == 1) || (m < M - 1 && m > 1 && board[m + 1][n] == 1 && board[m - 1][n] == 1 && board[m - 2][n] == 1) || (m > 0 && m < M - 2 && board[m + 2][n] == 1 && board[m + 1][n] == 1 && board[m - 1][n] == 1) || (m < M - 3 && board[m + 3][n] == 1 && board[m + 2][n] == 1 && board[m + 1][n] == 1))
									min = min - 20;
								if((board[m][n - 1] == 2 && board[m][n - 2] == 2 && board[m][n - 3] == 2 && n > 2) || (n > 1 && n < N - 1 && board[m][n + 1] == 2 && board[m][n - 1] == 2 && board[m][n - 2] == 2) || (n > 0 && n < N - 2 && board[m][n + 2] == 2 && board[m][n + 1] == 2 && board[m][n - 1] == 2) || (n < N - 3 && board[m][n + 3] == 2 && board[m][n + 2] == 2 && board[m][n + 1] == 2))
									max = max + 10;
								if((board[m][n - 1] == 1 && board[m][n - 2] == 1 && n > 2) || (n > 1 && n < N - 1 && board[m][n + 1] == 1 && board[m][n - 1] == 1) || (n > 0 && n < N - 2 && board[m][n + 2] == 1 && board[m][n + 1] == 1))
									min = min + 40;
								if((m > 2 && n < N - 3 && board[m - 1][n + 1] != 1 && board[m - 2][n + 2] != 1 && board[m - 3][n + 3] != 1) || (m > 1 && m < M - 1 && n > 0 && n < N - 2 && board[m + 1][n - 1] != 1 && board[m - 1][n + 1] != 1 && board[m - 2][n + 2] != 1) || (m > 0 && m < M - 2 && n > 1 && n < N - 1 && board[m + 2][n - 2] != 1 && board[m + 1][n - 1] != 1 && board[m - 1][n + 1] != 1) || (m < M - 3 && n > 2 && board[m + 3][n - 3] != 1 && board[m + 2][n - 2] != 1 && board[m + 1][n - 1] != 1))//���½ǵ����Ͻ�
									max++;
								if((m > 2 && n < N - 3 && board[m - 1][n + 1] != 2 && board[m - 2][n + 2] != 2 && board[m - 3][n + 3] != 2) || (m > 1 && m < M - 1 && n > 0 && n < N - 2 && board[m + 1][n - 1] != 2 && board[m - 1][n + 1] != 2 && board[m - 2][n + 2] != 2) || (m > 0 && m < M - 2 && n > 1 && n < N - 1 && board[m + 2][n - 2] != 2 && board[m + 1][n - 1] != 2 && board[m - 1][n + 1] != 2) || (m < M - 3 && n > 2 && board[m + 3][n - 3] != 2 && board[m + 2][n - 2] != 2 && board[m + 1][n - 1] != 2))//���½ǵ����Ͻ�
									min++;
								if((m > 2 && n > 2 && board[m - 1][n - 1] != 1 && board[m - 2][n - 2] != 1 && board[m - 3][n - 3] != 1) || (m > 1 && m < M - 1 && n > 1 && n < N - 1 && board[m + 1][n + 1] != 1 && board[m - 1][n - 1] != 1 && board[m - 2][n - 2] != 1) || (m > 0 && m < M - 2 && n > 0 && n < N - 2 && board[m + 2][n + 2] != 1 && board[m + 1][n + 1] != 1 && board[m - 1][n - 1] != 1) || (m < M - 3 && n < N - 3 && board[m + 3][n + 3] != 1 && board[m + 2][n + 2] != 1 && board[m + 1][n + 1] != 1))
									max++;
								if((m > 2 && n > 2 && board[m - 1][n - 1] != 2 && board[m - 2][n - 2] != 2 && board[m - 3][n - 3] != 2) || (m > 1 && m < M - 1 && n > 1 && n < N - 1 && board[m + 1][n + 1] != 2 && board[m - 1][n - 1] != 2 && board[m - 2][n - 2] != 2) || (m > 0 && m < M - 2 && n > 0 && n < N - 2 && board[m + 2][n + 2] != 2 && board[m + 1][n + 1] != 2 && board[m - 1][n - 1] != 2) || (m < M - 3 && n < N - 3 && board[m + 3][n + 3] != 2 && board[m + 2][n + 2] != 2 && board[m + 1][n + 1] != 2))
									min = min - 50;
							}
						}
					}*/
					//fp = max - min;
					
					//cout<<" getfp.....  "<<endl;
					//cout<<x1<<" "<<y1<<"****************"<<endl;
					fp = getfp(board, M, N, noX, noY);
					//cout<<"fp  ::::  "<<fp<<endl;
					if(fp < temp)
					{
						temp = fp;
						tempx = x1;
						tempy = y1;
					}
					if(temp < -10000000)
					    alphaB = false;
					//cout<<"temp :   "<<temp<<endl;
					board[x1][y1] = 0;
				}
			}
			//cout<<"����ֵ��  "<<temp<<"( "<<x <<" "<<y<<")"<<endl;
			if(temp > tempf)
			{
				tempf = temp;
				xx = x;
				yy = y;
			}
			//if(tempf > 100000000)
				//alphaB = false;
			board[x][y] = 0;
		}
	}
	x = xx;
	y = yy;//*****
	
//http://wenku.baidu.com/link?url=W4kELcy7AhucMx1bfhzPskRj3wg8rSetlxKhlQwGO3gyUnJeg04pzCjTfRuE5W-B_N7RQQSF1pWcXx1SIBLUplBxweupzYumAfwQpgtsl4S	
//alpha_beta��֦�㷨ʵ��******************
	if(alphaB){
    int *_top = new int[N];
	for(int t = 0; t < N; t++)
		_top[t] = top[t];
	//cout<<"cout top:::: "<<endl;
	//for(int t = 0; t < N; t++)
		//cout<<_top[t]<<" ";
	//cout<<"after cout top   "<<endl;
	node root;
	//cout<<"node *root:::::    "<<endl;
	root.father = NULL;
	//cout<<"father   "<<endl;
	root.isMax = true;
	//cout<<"isMax "<<endl;
	root.depth = depth;
	//cout<<"depth"<<endl;

	//root.x = lastX;
	//root.y = lastY;
	//cout<<"x   y"<<endl;
	//cout<<"b   ::::::"<<endl;
	for(int i = 0; i< M; i++)
	{
		for(int j = 0; j < N; j++)
			root.b[i][j] = board[i][j];
	}
	//cout<<"alpha   beta*************"<<endl;//OK!!!!!!
	node *_root;
	_root = &root;
	//cout<<_root->depth<<" "<<_root->isMax<<endl;//OK
    alphaBeta(_root, M, N, noX, noY, _top);//????????
	x = top[Y] - 1;
	y = Y;}
	//cout<<"x  :  "<<x<<"    y  :   "<<y<<endl;
//****************************************
	/*
		��Ҫ������δ���
	*/
	alphaB = true;
	clearArray(M, N, board);
	return new Point(x, y);
}


/*
	getPoint�������ص�Pointָ�����ڱ�dllģ���������ģ�Ϊ��������Ѵ���Ӧ���ⲿ���ñ�dll�е�
	�������ͷſռ䣬����Ӧ�����ⲿֱ��delete
*/
extern "C" __declspec(dllexport) void clearPoint(Point* p){
	delete p;
	return;
}

/*
	���top��board����
*/
void clearArray(int M, int N, int** board){
	for(int i = 0; i < M; i++){
		delete[] board[i];
	}
	delete[] board;
}


/*
	������Լ��ĸ�������������������Լ����ࡢ����������µ�.h .cpp�ļ�������ʵ������뷨
*/
