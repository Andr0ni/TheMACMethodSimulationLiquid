//---------------------------------------------------------------------------

#include <vcl.h>
#include <cmath>

#pragma hdrstop

#include "Unit2.h" 
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TForm2 *Form2;


using namespace std;

const double dt = 0.001;
double T_max = 60.0;
const double dx = 0.1;
const double dy = 0.1;
const int X = 2;
const int Y = 2;
const int cellX = 20;
const int cellY = 20;
const int a1 = 7;//y
const int a2 = 11;
const int b1 = 6;//x
const int b2 = 14;
const int num_marker = (a2 - a1 + 1) * (b2 - b1 + 1) * 4 - 1;
double Nu = 0.008;
double g[2] = {1, 0};//y,x

const int  Grid_size_x = 20;
const int Grid_size_y = 20;
const int cell_x  = 30;
const int cell_y = 30;

double grid[cellY + 2][cellX + 2];
double P[cellX + 2][cellY + 2];
double P0[cellX + 2][cellY + 2];
double U0[cellX + 2][cellY + 2];
double V0[cellX + 2][cellY + 2];
double U[cellX + 2][cellY + 2];
double V[cellX + 2][cellY + 2];
double D[cellX + 2][cellY + 2];
double R[cellX + 2][cellY + 2];
double Max = 0.;
double eps = 0.001;
double t = 0.;

bool stop=true;
TColor mColor;
TColor LiColor;

class Marker
{
public:
	Marker();
	~Marker();

	double x;
	double y;

	double xf, yf;

	double u, v;

	int i, j;

private:

};



Marker::Marker()
{
}

Marker::~Marker()
{
}

Marker M[num_marker + 1];

double fract(double x)
{
	return  abs(x-int(x));
}


double calc_u_i_j(int i, int j)
{
	return (U0[i + 1][j] + U0[i][j]) / 2.;
}

double calc_v_i_j(int i, int j)
{
	return (V0[i][j + 1] + V0[i][j]) / 2.;
}

double calc_v_u(int i, int j)
{
	return (U0[i][j] + U0[i][j - 1]) / 2. * (V0[i][j] + V0[i - 1][j]) / 2.;// äëÿ i-1/2,j-1/2
}

void start(double grid[cellY + 2][cellX + 2])
{
	int k = 0;
	for (int i = 0; i <= cellY + 1; i++)
	{
		for (int j = 0; j <= cellX + 1; j++)
		{
			P0[i][j] = 0.;
			U0[i][j] = 0.;
			V0[i][j] = 0.;
			P[i][j] = 0.;
			U[i][j] = 0.;
			V[i][j] = 0.;
			D[i][j] = 0.;
			R[i][j] = 0.;
        }
	}
	for (int i = 1; i < cellY; i++)
	{
		for (int j = 1; j < cellX; j++)
		{
			if (grid[i][j] == 2)
			{
				P0[i][j] = 1.;
				/*M[k].i = i - 2;//?????????????
				M[k].j = j - 2;
			    M[k].xf = 0.25 * dx;
				M[k].yf = 0.25 * dy;
				M[k].x = (j - 2) * dx + M[k].xf;
				M[k].y = (i - 2) * dy + M[k].yf;     */
				M[k].i = i - 1;//?????????????
				M[k].j = j - 1;
//				M[k].i = i;
//				M[k].j = j;
                M[k].x = 0.25*dx + M[k].j*dx;
                M[k].y = 0.25*dy + M[k].i*dy;
 //				M[k].xf = fract((M[k].x - int(M[k].x)-dx)/dx);
 //				M[k].yf = fract((M[k].y - int(M[k].y)-dy)/dy);
				M[k].xf = M[k].x / dx - M[k].j;
				M[k].yf = M[k].y / dy - M[k].i;

				M[k+1].i = i - 1;//?????????????
				M[k+1].j = j - 1;
//				M[k+1].i = i;
//				M[k+1].j = j;
                M[k+1].x = 0.25*dx + M[k+1].j*dx;
                M[k+1].y = 0.75*dy + M[k+1].i*dy;
//              M[k+1].xf = fract((M[k+1].x - int(M[k+1].x)-dx)/dx);
//				M[k+1].yf = fract((M[k+1].y - int(M[k+1].y)-dy)/dy);
				M[k+1].xf = M[k+1].x / dx - M[k+1].j;
				M[k+1].yf = M[k+1].y / dy - M[k+1].i;

				M[k+2].i = i - 1;//?????????????
				M[k+2].j = j - 1;
//				M[k+2].i = i;
//				M[k+2].j = j;
                M[k+2].x = 0.75*dx + M[k+2].j*dx;
                M[k+2].y = 0.25*dy + M[k+2].i*dy;
//              M[k+2].xf = fract((M[k+2].x - int(M[k+2].x)-dx)/dx);
//				M[k+2].yf = fract((M[k+2].y - int(M[k+2].y)-dy)/dy);
				M[k+2].xf = M[k+2].x / dx - M[k+2].j;
				M[k+2].yf = M[k+2].y / dy - M[k+2].i;

				M[k+3].i = i - 1;//?????????????
				M[k+3].j = j - 1;
//				M[k+3].i = i;
//				M[k+3].j = j;
                M[k+3].x = 0.75*dx + M[k+3].j*dx;
                M[k+3].y = 0.75*dy + M[k+3].i*dy;
//              M[k+3].xf = fract((M[k+3].x - int(M[k+3].x)-dx)/dx);
//				M[k+3].yf = fract((M[k+3].y - int(M[k+3].y)-dy)/dy);
				M[k+3].xf = M[k+3].x / dx - M[k+3].j;
				M[k+3].yf = M[k+3].y / dy - M[k+3].i;

			   /*	M[k + 1].i = i - 2;
				M[k + 1].j = j - 2;
				M[k + 1].xf = 0.25 * dx;
				M[k + 1].yf = 0.75 * dy;
				M[k + 1].x = (j - 2) * dx + M[k + 1].xf;
				M[k + 1].y = (i - 2) * dy + M[k + 1].yf;
				M[k + 2].i = i - 2;
				M[k + 2].j = j - 2;
				M[k + 2].xf = 0.75 * dx;
				M[k + 2].yf = 0.25 * dy;
				M[k + 2].x = (j - 2) * dx + M[k + 2].xf;
				M[k + 2].y = (i - 2) * dy + M[k + 2].yf;
				M[k + 3].i = i - 2;
				M[k + 3].j = j - 2;
				M[k + 3].xf = 0.75 * dx;
				M[k + 3].yf = 0.75 * dy;
				M[k + 3].x = (j - 2) * dx + M[k + 3].xf;
				M[k + 3].y = (i - 2) * dy + M[k + 3].yf;   */
				k += 4;
			}
		}
	}
}

//---------------------------------------------------------------------------
__fastcall TForm2::TForm2(TComponent* Owner)
	: TForm(Owner)
{
}
//---------------------------------------------------------------------------
void __fastcall TForm2::Button1Click(TObject *Sender)
{
if (stop)
{
   stop = false;
   Button1->Caption = "Stop";
}
else
{
   stop = true;
   Button1->Caption = "Continue";
}
while((t <= T_max) && (!stop))
	{
		t += dt;
		LabeledEdit2->Text = FloatToStrF(t,ffFixed,3,2);
		for (int i = 1; i < cellY; i++)
		{
			for (int j = 1; j < cellX; j++)
			{
				if (grid[i][j] == 2)
				{
					D[i][j] = (U0[i + 1][j] - U0[i][j]) / dx + (V0[i][j + 1] - V0[i][j]) / dy;
				}
			}
		}
		for (int i = 1; i < cellY; i++)
		{
			for (int j = 1; j < cellX; j++)
			{
				if (grid[i][j] == 2)
				{
					R[i][j] = (pow(calc_u_i_j(i + 1, j), 2) + pow(calc_u_i_j(i - 1, j), 2) - 2 * pow(calc_u_i_j(i, j), 2)) / dx / dx + (pow(calc_v_i_j(i, j + 1), 2) + pow(calc_v_i_j(i, j - 1), 2) - 2 * pow(calc_v_i_j(i, j), 2)) / dy / dy + 2 * (calc_v_u(i + 1, j + 1) + calc_v_u(i, j) - calc_v_u(i + 1, j) - calc_v_u(i, j + 1)) / dx / dy - D[i][j] / dt - Nu * ((D[i + 1][j] + D[i - 1][j] - 2 * D[i][j]) / dx / dx + (D[i][j + 1] + D[i][j - 1] - 2 * D[i][j]) / dy / dy);
				}
			}
		}
		do {
			Max = 0.;
			for (int i = 1; i < cellY; i++)
			{
				for (int j = 1; j < cellX; j++)
				{

					if (grid[i][j] == 2)
					{
						P[i][j] = ((P0[i + 1][j] + P0[i - 1][j]) / dy / dy + (P0[i][j + 1] + P0[i][j - 1]) / dx / dx + R[i][j]) / (2 / dx / dx + 2 / dy / dy);
						if (abs(P[i][j] - P0[i][j]) > Max)
						{
							Max = abs(P[i][j] - P0[i][j]);
						}
					}
					P[0][j]=P[1][j]+P[1][j]*g[0]*dy;
					P[cellY][j]=P[cellY - 1][j]-P[cellY - 1][j]*g[0]*dy;
					P[i][0]=P[i][1]-P[i][1]*g[1]*dx;
					P[i][cellX]=P[i][cellX - 1]+P[i][cellX - 1]*g[1]*dx;
				}
			}
			for (int i = 0; i <= cellY; i++)
			{
				for (int j = 0; j <= cellX; j++)
				{
					P0[i][j] = P[i][j];
				}
			}
		} while (Max > eps);

		for (int i = 0; i < cellY - 1; i++)//????? ò ê U[i + 1][j]
		{
			for (int j = 1; j < cellX; j++)
			{
				if (grid[i][j] == 2 || grid[i + 1][j] == 2 || grid[i - 1][j] == 2 || grid[i][j + 1] == 2 || grid[i][j - 1] == 2 /*|| grid[i + 1][j + 1] == 2 || grid[i - 1][j - 1] || grid[i - 1][j + 1] == 2 || grid[i + 1][j - 1] == 2*/)// ?????? íóæíî ëè óãëû ñ÷èòàòü
				{
					U[i + 1][j] = U0[i + 1][j] + dt * ((pow(calc_u_i_j(i, j), 2) - pow(calc_u_i_j(i + 1, j), 2)) / dx + (calc_v_u(i + 1, j) - calc_v_u(i + 1, j + 1)) / dy + g[0] + (P[i][j] - P[i + 1][j]) / dx + Nu * ((U0[i + 2][j] + U0[i][j] - 2 * U0[i + 1][j]) / dx / dx + (U0[i + 1][j + 1] + U0[i + 1][j - 1] - 2 * U0[i + 1][j]) / dy / dy));
				}
			}
		}

		for (int i = 1; i < cellY; i++)
		{
			for (int j = 0; j < cellX - 1; j++)
			{
				if (grid[i][j] == 2 || grid[i + 1][j] == 2 || grid[i - 1][j] == 2 || grid[i][j + 1] == 2 || grid[i][j - 1] == 2 /*|| grid[i + 1][j + 1] == 2 || grid[i - 1][j - 1] || grid[i - 1][j + 1] == 2 || grid[i + 1][j - 1] == 2*/)// ?????? íóæíî ëè óãëû ñ÷èòàòü
				{
					V[i][j + 1] = V0[i][j + 1] + dt * ((pow(calc_v_i_j(i, j), 2) - pow(calc_v_i_j(i, j + 1), 2)) / dy + (calc_v_u(i, j + 1) - calc_v_u(i + 1, j + 1)) / dx + g[1] + (P[i][j] - P[i][j + 1]) / dy + Nu * ((V0[i][j + 2] + V0[i][j] - 2 * V0[i][j + 1]) / dy / dy + (V0[i + 1][j + 1] + V0[i - 1][j + 1] - 2 * V0[i][j + 1]) / dx / dx));
				}
			}
		}

		for (int i = 1; i < cellY - 1; i++)
		{
			for (int j = 1; j < cellX; j++)
			{
				if (grid[i][j] == 0 && (grid[i + 1][j] == 2 || grid[i - 1][j] == 2))
				{
					if (grid[i - 1][j] == 0 && grid[i + 1][j] == 2)
					{
						U[i][j] = U[i + 1][j];
					}
					else if (grid[i + 1][j] == 0 && grid[i - 1][j] == 2)
					{
						U[i + 1][j] = U[i][j];
					}
				}
			}
		}
		for (int i = 1; i < cellY; i++)
		{
			for (int j = 1; j < cellX - 1; j++)
			{
				if (grid[i][j] == 0 && (grid[i][j + 1] == 2 || grid[i][j - 1] == 2))
				{
					if (grid[i][j - 1] == 0 && grid[i][j + 1] == 2)
					{
						V[i][j] = V[i][j + 1];
					}
					else if (grid[i][j + 1] == 0 && grid[i][j - 1] == 2)
					{
						V[i][j + 1] = V[i][j];
					}
				}
			}
		}

		/*for (int i = 1; i < cellY - 1; i++)
		{
			for (int j = 1; j < cellX; j++)
			{
				if (grid[i][j] == 0 && (grid[i + 1][j] == 2 || grid[i - 1][j] == 2 || grid[i][j + 1] == 2 || grid[i][j - 1] == 2))
				{
					if (grid[i + 1][j] == 2&&grid[i][j + 1] == 2&&grid[i][j - 1] == 2)
					{
						U[i][j] = U[i + 1][j];
					}
					else if (grid[i - 1][j] == 2&&grid[i][j + 1] == 2&&grid[i][j - 1] == 2)
					{
						U[i + 1][j] = U[i][j];
					}
				}
			}
		}
		for (int i = 1; i < cellY; i++)
		{
			for (int j = 1; j < cellX - 1; j++)
			{
				if (grid[i][j] == 0 && (grid[i + 1][j] == 2 || grid[i - 1][j] == 2 || grid[i][j + 1] == 2 || grid[i][j - 1] == 2))
				{
					if (grid[i + 1][j] == 2&&grid[i][j + 1] == 2&&grid[i - 1][j] == 2)
					{
						V[i][j] = V[i][j + 1];
					}
					else if (grid[i + 1][j] == 2&&grid[i][j - 1] == 2&&grid[i - 1][j] == 2)
					{
						V[i][j + 1] = V[i][j];
					}
				}
			}
		}*/

		for(int j=0;j<=cellX;j++)
		{
			V[0][1]=V[1][1];
			V[cellY][0]=V[cellY - 1][0];
			U[cellY][j]= 0;
			U[cellY + 1][j]= -U[cellY - 1][j];

			V[cellY][j + 1]=V[cellY - 1][j + 1];
			U[1][j]= 0;
			U[0][j]= -U[2][j];

			V[0][j + 1]=V[1][j + 1];
		}
		for(int i=0;i<=cellY;i++)
		{
			U[1][0]=U[1][1];
			U[1][cellX]=U[1][cellX - 1];
			V[i][cellX]= 0;
			V[i][cellX + 1]= -V[i][cellX - 1];

			U[i + 1][0]=U[i + 1][1];
			V[i][1]= 0;
			V[i][0]= -V[i][2];

			U[i + 1][cellX]=U[i + 1][cellX - 1];
		}

		for (int k = 0; k <= num_marker; k++)
		{
//			if (M[k].xf > 0.5)//åñëè ìàðêåð íàõîäèòñÿ ïðàâåå
//			{
//				M[k].x = M[k].x + (V[M[k].i][M[k].j] + V[M[k].i][M[k].j+1] /* + V[M[k].i + 1][M[k].j] + V[M[k].i + 1][M[k].j + 1]*/) * dt / 2.;
//			}
//			else//ëåâåå
//			{
//				M[k].x = M[k].x + (V[M[k].i][M[k].j] + V[M[k].i][M[k].j-1] /* + V[M[k].i - 1][M[k].j] + V[M[k].i - 1][M[k].j + 1]*/) * dt / 2.;
//			}
//			if (M[k].yf > 0.5)//ââåðõó
//			{
//				M[k].y = M[k].y + (U[M[k].i][M[k].j] + U[M[k].i-1 ][M[k].j] /* + U[M[k].i][M[k].j + 1] + U[M[k].i + 1][M[k].j + 1]*/) * dt / 2.;
//			}
//			else//âíèçó
//			{
//				M[k].y = M[k].y + (U[M[k].i][M[k].j] + U[M[k].i+1][M[k].j ] /* + U[M[k].i][M[k].j + 1] + U[M[k].i + 1][M[k].j + 1]*/) * dt / 2.;
//			}
			if (M[k].yf >= 0.5)// åñëè ìàðêåð íàõîäèòñÿ  ñíèçó
			{
				M[k].x = M[k].x + (V[M[k].i][M[k].j] + V[M[k].i][M[k].j + 1] + V[M[k].i + 1][M[k].j] + V[M[k].i + 1][M[k].j + 1]) * dt / 4.;
			}
			else// ñâåðõó
			{
				M[k].x = M[k].x + (V[M[k].i][M[k].j] + V[M[k].i][M[k].j + 1] + V[M[k].i - 1][M[k].j] + V[M[k].i - 1][M[k].j + 1]) * dt / 4.;
			}
			if (M[k].xf >= 0.5)// ïðàâåå
			{
				M[k].y = M[k].y + (U[M[k].i][M[k].j] + U[M[k].i + 1][M[k].j] + U[M[k].i][M[k].j + 1] + U[M[k].i + 1][M[k].j + 1]) * dt / 4.;
			}
			else// ëåâåå
			{
				M[k].y = M[k].y + (U[M[k].i][M[k].j] + U[M[k].i + 1][M[k].j] + U[M[k].i][M[k].j - 1] + U[M[k].i + 1][M[k].j - 1]) * dt / 4.;
			}
			M[k].i = int(M[k].y / dy + 1);
			M[k].j = int(M[k].x / dx + 1);
			M[k].xf = M[k].x / dx + 1 - M[k].j;
			M[k].yf = M[k].y / dy + 1 - M[k].i;
		}

		for (int i = 0; i <= cellY; i++)
		{
			for (int j = 0; j <= cellX; j++)
			{
				U0[i][j] = U[i][j];
				V0[i][j] = V[i][j];
			}
		}
        for (int i = 0; i <= cellY; i++)
		{
			for (int j = 0; j <= cellX; j++)
			{
				D[i][j] = 0.;
				R[i][j] = 0.;
				P[i][j] = 0.;
				//U[i][j] = 0.;  //???????âðîäå êàê íóæíî çàíóëÿòü
				//V[i][j] = 0.;
			}
		}
		for (int i = 1; i < cellY; i++)
		{
			for (int j = 1; j < cellX; j++)
			{
				grid[i][j] = 0;
			}
		}
		for (int k = 0; k <= num_marker; k++)
		{
			grid[M[k].i][M[k].j] = 2;
		}

	for (int i = 1; i < cellY; i++)
	{
		for (int j = 1; j < cellX; j++)
		{
            if (grid[i][j] == 2){
             if(CheckBox1->Checked){
                  Image1->Canvas->Brush->Color = LiColor;
             }
             else {
                  Image1->Canvas->Brush->Color = clWhite;
             }
            Image1->Canvas->Pen->Color = LiColor;
            }
			else{
				 if (grid[i][j] == -1){
                 Image1->Canvas->Brush->Color = clGray;
                 Image1->Canvas->Pen->Color = clGray;
                 }
				 else {
				 Image1->Canvas->Brush->Color = clWhite;
				 Image1->Canvas->Pen->Color = clGray;
				 }
			}
			Image1->Canvas->Rectangle(cell_x*j,cell_y*i,cell_x*(j+1),cell_y*(i+1));
		}
	}
	int r=2;
	for (int k = 0; k <= num_marker; k++){
		if(M[k].i!=0&&M[k].i!=cellY&&M[k].j!=0&&M[k].j!=cellX)
		{
			Image1->Canvas->Brush->Color = mColor;
			Image1->Canvas->Pen->Color = mColor;
			//Image1->Canvas->Pixels[M[k].xf*cell_x+M[k].i*cell_x][(M[k].yf)*cell_y+M[k].j*cell_y] = clLime;
			Image1->Canvas->Ellipse(M[k].xf*cell_x+M[k].j*cell_x-r,M[k].yf*cell_y+M[k].i*cell_y-r,M[k].xf*cell_x+M[k].j*cell_x+r,M[k].yf*cell_y+M[k].i*cell_y+r);
		//    Image1->Canvas->Pixels[M[k].xf*cell_x+M[k].i*cell_x][M[k].yf*cell_y+M[k].j*cell_y] = clLime;
		}
	}
	for (int i = 0; i <= cellY + 1; i++)
	{
		for (int j = 0; j <= cellX + 1; j++)
		{
			if (i == 0 || i == cellY || j == 0 || j == cellX)
			{
                Image1->Canvas->Brush->Color = clGray;
				Image1->Canvas->Pen->Color = clGray;
			}
		}
	}
       Application->ProcessMessages();
	}

}
//---------------------------------------------------------------------------
//---------------------------------------------------------------------------


void __fastcall TForm2::Button2Click(TObject *Sender)
{
t=0;
LabeledEdit2->Text = FloatToStrF(t,ffFixed,3,2);
Button1->Enabled = true;
mColor = ColorBox1->Selected;
LiColor = ColorBox2->Selected;
T_max = StrToFloat(LabeledEdit1->Text);
g[0] = StrToFloat(LabeledEdit3->Text);
g[1] = StrToFloat(LabeledEdit4->Text);
Nu = StrToFloat(LabeledEdit5->Text);



	for (int i = 0; i <= cellY; i++)//«àäàåì ñåòêó
	{
		for (int j = 0; j <= cellX; j++)
		{
			if (i == 0 || i == cellY || j == 0 || j == cellX)
			{
				grid[i][j] = -1;
			}
			else
			{
				grid[i][j] = 0;
			}
		}
	}
	for (int i = a1; i <= a2; i++) // âîäà
	{
		for (int j = b1; j <= b2; j++)
		{
			grid[i][j] = 2;
		}
	}
    /*
	for (int i = 0; i <= cellY; i++)
	{
		for (int j = 0; j <= cellX; j++)
		{
			cout << grid[i][j] << "\t";
		}
		//cout << endl;
	}     */


	Image1->Height = cell_y*(cellY+1);
    Image1->Width =  cell_x*(cellX+1);

    //Form2->Width =  cell_x*(cellX+1)+Button1->Width*3;

	for (int i = 0; i <= cellY; i++)
	{
		for (int j = 0; j <= cellX; j++)
		{
            if (grid[i][j] == 2){
             if(CheckBox1->Checked){
                  Image1->Canvas->Brush->Color = LiColor;
             }
             else {
                  Image1->Canvas->Brush->Color = clWhite;
             }
            Image1->Canvas->Pen->Color = LiColor;
            }
            else{
				 if (grid[i][j] == -1){
				 Image1->Canvas->Brush->Color = clGray;
				 Image1->Canvas->Pen->Color = clGray;
                 }
				 else {
                 Image1->Canvas->Brush->Color = clWhite;
                 Image1->Canvas->Pen->Color = clGray;
				 }
            }
			Image1->Canvas->Rectangle(cell_x*j,cell_y*i,cell_x*(j+1),cell_y*(i+1));
		}
	}
	start(grid);//íà÷àëüíûå óñëîâè¤

	int r=2;
	for (int k = 0; k <= num_marker; k++){
    Image1->Canvas->Brush->Color = mColor;
    Image1->Canvas->Pen->Color = mColor;
	Image1->Canvas->Ellipse(M[k].xf*cell_x+(M[k].j+1)*cell_x-r,M[k].yf*cell_y+(M[k].i+1)*cell_y-r,M[k].xf*cell_x+(M[k].j+1)*cell_x+r,M[k].yf*cell_y+(M[k].i+1)*cell_y+r);
    }
}
//---------------------------------------------------------------------------

void __fastcall TForm2::ColorBox1Change(TObject *Sender)
{
  mColor = ColorBox1->Selected;
}
//---------------------------------------------------------------------------

