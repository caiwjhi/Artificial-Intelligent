#include <iostream>
#include <memory>
#include "Point.h"
#include "Strategy.h"

/*#define dead1 0;
#define alive1 2;
#define dead2 3;
#define alive2 100;
#define dead3 500000;
#define alive3 500000;//棋局的各种状态的权值，死1：只有一个棋子连起来，旁边还有一个对方棋子堵着；活2：两个棋子连起来，两端没有对方棋子堵着。。
*/
const int dead1 = 10;
const int alive1 = 30;
const int dead2 = 50;
const int alive2 = 8000;
const int dead3 = 9000;
const int alive3 = 100000;
const int dead4 = 20000000;//死活都一样
const int alive4 = 20000000;
const int depth = 5;//alpha_beta算法搜寻树的高度，深度
//const int legth = 50000;//(int)((pow(12.0, 5*1.0 + 1) - 1) / 11);//节点最多数目
using namespace std;

int Y = -1;
bool alphaB = true;
struct node
{
public:
    node *father;
	//node *sons[12];//最多有12个孩子
	int alpha;//极大点的下界
	int beta;//极小点的上界
	bool isMax;//极大或极小
   // bool hasfp;//权值
	int depth;//深度
	int **b;//该点的棋盘状态
	//int x, y;//到达该点所改变的棋子坐标
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
		//hasfp = false;//叶节点计算fp。其他节点只须计算上下界，也可用来判断是否已计算过上下界，也可以判断现在是处于上溯还是下溯
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



int *getND(int **b, int max, int M, int N, int noX, int noY)//M N 是棋局长度和宽度
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
	if(max == 1)//计算我方的状态数
	{
		//cout<<"max == 1  "<<endl;
		for(i = 0; i < M; i++)//横着，一行一行计算
		{
			for(j = 0; j < N; j++)
			{
				if(b[i][j] == 2)//遇到我方棋子
				{
					if(j == N - 1)
						num[0]++;//死一
					else
					{
						if(j == 0)//死
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
							if(b[i][j - 1] == 1 || (i == noX && j - 1 == noY))//对方棋子或该处不可下子堵着，为死,,,,此处考虑了不可落子点
							{
								num[2*(k - j) - 2]++;
							}
							else
							{
								if((k < N) && (b[i][k] == 0) && !(i == noX && k == noY))//右边界也空着，活
									num[2*(k - j) - 1]++;
								else
									num[2*(k - j) - 2]++;//死
							}
						}
					}
				}
			}
		}
		//cout<<"after 一行一行：：："<<endl;
		/*for(int t = 0; t < 10; t++)
		{
			cout<<num[t]<<endl;
		}*/
		for(j = 0; j < N; j++)
		{
			for(i = 0; i < M; i++)//一列一列
			{
				if(b[i][j] == 2)//遇到我方棋子
				{
					if(i == M - 1)
						num[0]++;//死一
					else
					{
						if(i == 0)//死
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
							if(b[i - 1][j] == 1 || (i - 1 == noX && j == noY))//对方棋子或该处不可下子堵着，为死,,,,此处考虑了不可落子点
							{
								num[2*(k - i) - 2]++;
							}
							else
							{
								if((k < M) && (b[k][j] == 0) && !(k == noX && j == noY))//右边界也空着，活
									num[2*(k - i) - 1]++;
								else
									num[2*(k - i) - 2]++;//死
							}
						}
					}
				}
			}
		}
		for(i = 0; i < M; i++)//斜着， 第一次找到的是一串里的第一个，死一活一均不计
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
						for(k = 0; k + j < N && k + i < M; k++)//从左上角到右下角
						{
							if(b[i + k][j + k] != 2)
							{
								break;
							}
							if((i > 0 && j > 0) && b[i - 1][j - 1] == 0 && !(i - 1 == noX && j - 1 == noY) && i + k < M && j + k < N && b[i + k][j + k] == 0 && !(i + k == noX && j + k == noY))//活的条件
							{
								//if(k > 1)
									num[2 * k - 1]++;
							}
							else//死
							{
								//if(k > 1)//去掉单独一个棋子的情况
									num[2*k - 2]++;
							}
						}
						for(k = 0; i + k < M && j - k >= 0; k++)
						{
							if(b[i + k][j - k] != 2)
								break;
							if(i > 0 && j < N - 1 && b[i - 1][j + 1] == 0 && !(i - 1 == noX && j + 1 == noY) && i + k < M && j - k >= 0 && b[i + k][j - k] == 0 && !(i + k == noX && j - k == noY))//活的条件
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
		for(i = 0; i < M; i++)//横着，一行一行计算
		{
			for(j = 0; j < N; j++)
			{
				if(b[i][j] == 1)//遇到我方棋子
				{
					if(j == N - 1)
						num[0]++;//死一
					else
					{
						if(j == 0)//死
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
							if(b[i][j - 1] == 2 || (i == noX && j - 1 == noY))//对方棋子或该处不可下子堵着，为死,,,,此处考虑了不可落子点
							{
								num[2*(k - j) - 2]++;
							}
							else
							{
								if((k < N) && (b[i][k] == 0) && !(i == noX && k == noY))//右边界也空着，活
									num[2*(k - j) - 1]++;
								else
									num[2*(k - j) - 2]++;//死
							}
						}
					}
				}
			}
		}
		for(j = 0; j < N; j++)
		{
			for(i = 0; i < M; i++)//一列一列
			{
				if(b[i][j] == 1)//遇到我方棋子
				{
					if(i == M - 1)
						num[0]++;//死一
					else
					{
						if(i == 0)//死
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
							if(b[i - 1][j] == 2 || (i - 1 == noX && j == noY))//对方棋子或该处不可下子堵着，为死,,,,此处考虑了不可落子点
							{
								num[2*(k - i) - 2]++;
							}
							else
							{
								if((k < M) && (b[k][j] == 0) && !(k == noX && j == noY))//右边界也空着，活
									num[2*(k - i) - 1]++;
								else
									num[2*(k - i) - 2]++;//死
							}
						}
					}
				}
			}
		}
		for(i = 0; i < M; i++)//斜着， 第一次找到的是一串里的第一个，死一活一均不计
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
						for(k = 0; k + j < N && k + i < M; k++)//从左上角到右下角
						{
							if(b[i + k][j + k] != 1)
							{
								break;
							}
							if((i > 0 && j > 0) && b[i - 1][j - 1] == 0 && !(i - 1 == noX && j - 1 == noY) && i + k < M && j + k < N && b[i + k][j + k] == 0 && !(i + k == noX && j + k == noY))//活的条件
							{
								//if(k > 1)
									num[2 * k - 1]++;
							}
							else//死
							{
								//if(k > 1)//去掉单独一个棋子的情况
									num[2*k - 2]++;
							}
						}
						for(k = 0; i + k < M && j - k >= 0; k++)
						{
							if(b[i + k][j - k] != 1)
								break;
							if(i > 0 && j < N - 1 && b[i - 1][j + 1] == 0 && !(i - 1 == noX && j + 1 == noY) && i + k < M && j - k >= 0 && b[i + k][j - k] == 0 && !(i + k == noX && j - k == noY))//活的条件
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

int getfp(int **b, int M, int N, int noX, int noY)//该函数根据当前棋局状态计算该状态下的评估值fp = max - min
{
	/*for(int m = 0; m < M; m++)
	{
		for(int n = 0; n < N; n++)
		{
			cout<<b[m][n]<<" ";
		}
		cout<<endl;
	}*///OK!!!!
	int *numOfADmax;//各种状态的数量，状态死1下标是0，状态活1下标是1，。。。
	numOfADmax = new int[8];
	memset(numOfADmax, 0, sizeof(numOfADmax));
	//cout<<"getND   "<<endl;
	numOfADmax = getND(b, 1, M, N, noX, noY);//我/程序/max 方的状态数量
	//cout<<"我方的状态数量已知"<<endl;
	//for(int t = 0; t < 8; t++)
		//cout<<numOfADmax[t]<<" ";
	//cout<<"after cout numOfADmax "<<endl;
	int *numOfADmin;
	numOfADmin = new int[8];
	numOfADmin = getND(b, 0, M, N, noX, noY);//对方的状态数量
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
	if(root->depth == 0)//叶节点
		return getfp(root->b, M, N, noX, noY);
	int i, j;
	int x, y;
	int value;//需要被返回的值，
	//cout<<"root :    "<<endl;
	//cout<<root->depth<<"  "<<root->isMax<<endl;
	for(j = 0; j < N; j++)
	{
		if((_top[j] - 1 >= 0) && (_top[j] - 1 < M) && !(_top[j] - 1 == noX && j == noY))//该列可以下子
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
			next->depth = root->depth - 1;//深度在减少
			//cout<<"root ....."<<endl;
			if(root->isMax)//父节点为极大节点
			{
				//cout<<"isMax  ...."<<endl;
				next->isMax = false;

				next->b[x][y] = 2;//极大节点下子变成极小节点
				if(j == 0)//最左边的孩子
				{
					/*if(getfp(next->b, M, N, noX, noY) > 1000000000 && root->depth == depth)//必胜
					{	
						Y = j;
						cout<<"getfp(next->b, M, N, noX, noY)   "<<getfp(next->b, M, N, noX, noY)<<endl;
						break;
					}*/
					//cout<<"最左边.....  "<<endl;
					next->beta = alphaBeta(next, M, N, noX, noY, _top);
					//cout<<"递归之后，next->beta    "<<next->beta<<endl;
					root->alpha = next->beta;
					if(root->depth == depth)
						Y = j;
					//cout<<"Y    ::   "<<Y<<endl;
					value = next->beta;
				}
				else
				{
					//cout<<"next is not the first child  ......."<<endl;
					/*if(getfp(next->b, M, N, noX, noY) > 1000000000 && root->depth == depth)//必胜
					{	
						Y = j;
						cout<<"必胜  ..."<<getfp(next->b, M, N, noX, noY)<<endl;
						break;
					}*/
					//node _temp;
					node *temp = root->father;
					//temp = root->father;
					while(temp)
					{
						if(root->alpha >= temp->beta)//beta剪枝
						{
							break;//遍历停止
						}
						if(temp->father)
						   temp = temp->father->father;
						else  temp=NULL;
					}
					if(temp)//剪枝了
					{
						//cout<<"剪枝   。。。。。"<<endl;
						value = root->alpha;
						//continue;
						_top[j]++;
						break;
					}
					else//????没有返回
					{
						//cout<<"没有剪枝"<<endl;
						next->beta = alphaBeta(next, M, N, noX, noY, _top);//向下拓展
						//cout<<"next -> beta  .....  "<<next->beta<<endl;
						//root->father->beta = root->alpha;//向上更新 
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
				/*if(getfp(next->b, M, N, noX, noY) < -1000000 && root->depth == depth - 1)//必输
					{	
						value = getfp(next->b, M, N, noX, noY);
						cout<<"输了。。。。"<<endl;
						break;
					}*/
				if(j == 0)//最左边的孩子
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
						if(root->beta <= temp->alpha)//beta剪枝
						{
							break;//遍历停止
						}
						if(temp->father)
						  temp = temp->father->father;
						else
							temp = NULL;
					}
					if(temp)//剪枝了
					{
						//cout<<"剪枝  。。。。。  "<<endl;
						value = root->beta;
						//continue;
						_top[j]++;
						break;

					}
					else//??????
					{
						next->alpha = alphaBeta(next, M, N, noX, noY, _top);//向下拓展
						//cout<<"next alpha  ....  "<<next->alpha<<endl;
						/*if(root->depth == depth)
						{
							//root->father->alpha = root->beta;//向上更新
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
		//还原
		}
	}
	return value;
}
/*
	策略函数接口,该函数被对抗平台调用,每次传入当前状态,要求输出你的落子点,该落子点必须是一个符合游戏规则的落子点,不然对抗平台会直接认为你的程序有误
	
	input:
		为了防止对对抗平台维护的数据造成更改，所有传入的参数均为const属性
		M, N : 棋盘大小 M - 行数 N - 列数 均从0开始计， 左上角为坐标原点，行用x标记，列用y标记
		top : 当前棋盘每一列列顶的实际位置. e.g. 第i列为空,则_top[i] == M, 第i列已满,则_top[i] == 0//如果下次在此列下子，则下在_top[i] - 1处
		_board : 棋盘的一维数组表示, 为了方便使用，在该函数刚开始处，我们已经将其转化为了二维数组board
				你只需直接使用board即可，左上角为坐标原点，数组从[0][0]开始计(不是[1][1])
				board[x][y]表示第x行、第y列的点(从0开始计)
				board[x][y] == 0/1/2 分别对应(x,y)处 无落子/有用户的子/有程序的子,不可落子点处的值也为0
		lastX, lastY : 对方上一次落子的位置, 你可能不需要该参数，也可能需要的不仅仅是对方一步的
				落子位置，这时你可以在自己的程序中记录对方连续多步的落子位置，这完全取决于你自己的策略
		noX, noY : 棋盘上的不可落子点(注:其实这里给出的top已经替你处理了不可落子点，也就是说如果某一步
				所落的子的上面恰是不可落子点，那么UI工程中的代码就已经将该列的top值又进行了一次减一操作，
				所以在你的代码中也可以根本不使用noX和noY这两个参数，完全认为top数组就是当前每列的顶部即可,
				当然如果你想使用lastX,lastY参数，有可能就要同时考虑noX和noY了)
		以上参数实际上包含了当前状态(M N _top _board)以及历史信息(lastX lastY),你要做的就是在这些信息下给出尽可能明智的落子点
	output:
		你的落子点Point
*/


extern "C" __declspec(dllexport) Point* getPoint(const int M, const int N, const int* top, const int* _board, 
	const int lastX, const int lastY, const int noX, const int noY){
	/*
		不要更改这段代码
	*/
	int x = -1, y = -1;//最终将你的落子点存到x,y中
	int** board = new int*[M];
	for(int i = 0; i < M; i++){
		board[i] = new int[N];
		for(int j = 0; j < N; j++){
			board[i][j] = _board[i * N + j];
		}
	}
	
	/*
		根据你自己的策略来返回落子点,也就是根据你的策略完成对x,y的赋值
		该部分对参数使用没有限制，为了方便实现，你可以定义自己新的类、.h文件、.cpp文件
	*/
	//Add your own code below
	/*
     //a naive example
	for (int i = N-1; i >= 0; i--) {
		if (top[i] > 0) {
			x = top[i] - 1;
			y = i;
			break;
		}//直接从右边一次下子
	}*/
    //极大极小*****
	
	int i, j;
	int tempf = -1000000000;//最小，求最大
	int xx, yy;//记录最大tempf时的坐标x, y
	for(i = 0; i < N; i++)//列
	{
		if(top[i] > 0)
		{
			x = top[i] - 1;
			y = i;
			board[x][y] = 2;
			int temp = 2000000000;//初始化最大，求最小
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
				if(y1 > -1)//在该列找到可落子点
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
								if((m > 2 && board[m - 1][n] != 1 && board[m - 2][n] != 1 && board[m - 3][n] != 1) || (m < M - 1 && m > 1 && board[m + 1][n] != 1 && board[m - 1][n] != 1 && board[m - 2][n] != 1) || (m > 0 && m < M - 2 && board[m + 2][n] != 1 && board[m + 1][n] != 1 && board[m - 1][n] != 1) || (m < M - 3 && board[m + 3][n] != 1 && board[m + 2][n] != 1 && board[m + 1][n] != 1))//从下到上
									max++;
								if((m > 2 && board[m - 1][n] != 2 && board[m - 2][n] != 2 && board[m - 3][n] != 2) || (m > 1 && m < M - 1 && board[m + 1][n] != 2 && board[m - 1][n] != 2 && board[m - 2][n] != 2) || (m > 0 && m < M - 2 && board[m + 2][n] != 2 && board[m + 1][n] != 2 && board[m - 1][n] != 2) || (m < M - 3 && board[m + 3][n] != 2 && board[m + 2][n] != 2 && board[m + 1][n] != 2))//从下到上
									min--;
								if((board[m][n - 1] != 1 && board[m][n - 2] != 1 && board[m][n - 3] != 1 && n > 2) || (n > 1 && n < N - 1 && board[m][n + 1] != 1 && board[m][n - 1] != 1 && board[m][n - 2] != 1) || (n > 0 && n < N - 2 && board[m][n + 2] != 1 && board[m][n + 1] != 1 && board[m][n - 1] != 1) || (n < N - 3 && board[m][n + 3] != 1 && board[m][n + 2] != 1 && board[m][n + 1] != 1))//从左到右
									max++;
								if((board[m][n - 1] != 2 && board[m][n - 2] != 2 && board[m][n - 3] != 2 && n > 2) || (n > 1 && n < N - 1 && board[m][n + 1] != 2 && board[m][n - 1] != 2 && board[m][n - 2] != 2) || (n > 0 && n < N - 2 && board[m][n + 2] != 2 && board[m][n + 1] != 2 && board[m][n - 1] != 2) || (n < N - 3 && board[m][n + 3] != 2 && board[m][n + 2] != 2 && board[m][n + 1] != 2))//从左到右，横着
									min--;
								if((m > 2 && board[m - 1][n] == 2 && board[m - 2][n] == 2 && board[m - 3][n] == 2) || (m < M - 1 && m > 1 && board[m + 1][n] == 2 && board[m - 1][n] == 2 && board[m - 2][n] == 2) || (m > 0 && m < M - 2 && board[m + 2][n] == 2 && board[m + 1][n] == 2 && board[m - 1][n] == 2) || (m < M - 3 && board[m + 3][n] == 2 && board[m + 2][n] == 2 && board[m + 1][n] == 2))
									max = max + 10;
								if((m > 2 && board[m - 1][n] == 1 && board[m - 2][n] == 1 && board[m - 3][n] == 1) || (m < M - 1 && m > 1 && board[m + 1][n] == 1 && board[m - 1][n] == 1 && board[m - 2][n] == 1) || (m > 0 && m < M - 2 && board[m + 2][n] == 1 && board[m + 1][n] == 1 && board[m - 1][n] == 1) || (m < M - 3 && board[m + 3][n] == 1 && board[m + 2][n] == 1 && board[m + 1][n] == 1))
									min = min - 20;
								if((board[m][n - 1] == 2 && board[m][n - 2] == 2 && board[m][n - 3] == 2 && n > 2) || (n > 1 && n < N - 1 && board[m][n + 1] == 2 && board[m][n - 1] == 2 && board[m][n - 2] == 2) || (n > 0 && n < N - 2 && board[m][n + 2] == 2 && board[m][n + 1] == 2 && board[m][n - 1] == 2) || (n < N - 3 && board[m][n + 3] == 2 && board[m][n + 2] == 2 && board[m][n + 1] == 2))
									max = max + 10;
								if((board[m][n - 1] == 1 && board[m][n - 2] == 1 && n > 2) || (n > 1 && n < N - 1 && board[m][n + 1] == 1 && board[m][n - 1] == 1) || (n > 0 && n < N - 2 && board[m][n + 2] == 1 && board[m][n + 1] == 1))
									min = min + 40;
								if((m > 2 && n < N - 3 && board[m - 1][n + 1] != 1 && board[m - 2][n + 2] != 1 && board[m - 3][n + 3] != 1) || (m > 1 && m < M - 1 && n > 0 && n < N - 2 && board[m + 1][n - 1] != 1 && board[m - 1][n + 1] != 1 && board[m - 2][n + 2] != 1) || (m > 0 && m < M - 2 && n > 1 && n < N - 1 && board[m + 2][n - 2] != 1 && board[m + 1][n - 1] != 1 && board[m - 1][n + 1] != 1) || (m < M - 3 && n > 2 && board[m + 3][n - 3] != 1 && board[m + 2][n - 2] != 1 && board[m + 1][n - 1] != 1))//左下角到右上角
									max++;
								if((m > 2 && n < N - 3 && board[m - 1][n + 1] != 2 && board[m - 2][n + 2] != 2 && board[m - 3][n + 3] != 2) || (m > 1 && m < M - 1 && n > 0 && n < N - 2 && board[m + 1][n - 1] != 2 && board[m - 1][n + 1] != 2 && board[m - 2][n + 2] != 2) || (m > 0 && m < M - 2 && n > 1 && n < N - 1 && board[m + 2][n - 2] != 2 && board[m + 1][n - 1] != 2 && board[m - 1][n + 1] != 2) || (m < M - 3 && n > 2 && board[m + 3][n - 3] != 2 && board[m + 2][n - 2] != 2 && board[m + 1][n - 1] != 2))//左下角到右上角
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
			//cout<<"评估值：  "<<temp<<"( "<<x <<" "<<y<<")"<<endl;
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
//alpha_beta剪枝算法实现******************
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
		不要更改这段代码
	*/
	alphaB = true;
	clearArray(M, N, board);
	return new Point(x, y);
}


/*
	getPoint函数返回的Point指针是在本dll模块中声明的，为避免产生堆错误，应在外部调用本dll中的
	函数来释放空间，而不应该在外部直接delete
*/
extern "C" __declspec(dllexport) void clearPoint(Point* p){
	delete p;
	return;
}

/*
	清除top和board数组
*/
void clearArray(int M, int N, int** board){
	for(int i = 0; i < M; i++){
		delete[] board[i];
	}
	delete[] board;
}


/*
	添加你自己的辅助函数，你可以声明自己的类、函数，添加新的.h .cpp文件来辅助实现你的想法
*/
