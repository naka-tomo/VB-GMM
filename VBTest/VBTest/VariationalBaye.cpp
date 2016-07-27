#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include "VariationalBayes.h"
#include <float.h>
#include <windows.h>
#include "gamma.h"
#include "utility.h"

#pragma warning(disable: 4996)

// 0～1の乱数を返す
double RandF()
{
	unsigned int rnd;
	rand_s( &rnd );
	return (double)rnd/UINT_MAX;
}

// 整数の乱数を返す
unsigned int RandD()
{
	unsigned int rnd;
	rand_s( &rnd );
	return rnd;
}

CVariationalBayes::CVariationalBayes()
{
	//構造体初期化
	memset( m_models , 0 , MAXMODEL*sizeof(MODEL) );

	for(int i=0;i<DIM;i++)
	{
		m_data[i] = NULL;
	}
}


CVariationalBayes::~CVariationalBayes()
{
	FreeMemory();
}



// プログラム全体を通じて最後に一度実行される関数
void CVariationalBayes::FreeMemory()
{
	for(int i=0 ; i<MAXMODEL ; i++) FreeModel( &(m_models[i]) );
}

// モデル構造体を確保
int CVariationalBayes::AllocModel(MODEL *model , int modelnum )
{
	model->m = new PARAMETER[modelnum+1];

	for(int j=0 ; j<=modelnum ;j++)
	{
		model->m[j].z = new double[m_numData];
	}

	model->m_freeEnergies = 0;
	model->ModelNum = modelnum+1;

	return 0;
}

int CVariationalBayes::FreeModel( MODEL *model )
{
	if( model->m )
	{
		for(int j=0 ; j<model->ModelNum ;j++)
		{
			delete [] model->m[j].z;
			model->m[j].z = NULL;
		}

		delete [] model->m;
		model->m = NULL;
	}

	model->m_freeEnergies = 0;
	model->ModelNum = 0;

	return 0;
}

// 学習関数
void CVariationalBayes::Learn( std::string SaveDir )
{
	//データ保存用ディレクトリ
	if( SaveDir.size() )
	{
		if( *(SaveDir.end()-1) != '/' || *(SaveDir.end()-1) != '\\' ) SaveDir += "/";
	}
	else
	{
		printf("ディレクトリを指定してください。");
		return;
	}
	CreateDirectory( SaveDir.c_str() , NULL );


	//メモリ確保
	for(int i=0 ; i<MAXMODEL ; i++)
	{
		if( m_models[i].m ) FreeModel( &(m_models[i]) );
		AllocModel( &(m_models[i]) , i );
		m_freeEnergies[i] = -DBL_MAX;
	}


	//model数を変えて最適事後分布を計算
	for(int M=1;M<MAXMODEL;M++)
	{
		printf("モデル数 %02d 計算開始" , M+1 );

		PARAMETER *params = m_models[M].m;

		for(int it=0 ; it<ITERATE ; it++ )
		{
			double oldF = -DBL_MAX;
			double Fm = 0;
			
			//Step1-1:初期パラメータ決定
			InitializeParameter(M);

			//Step1-2:変数初期値決定
			for(int m=0;m<=M;m++){
				SetInitializingParameter(&params[m],M,m);
			}

			for(int loop=0 ; ; loop++ )
			{
				//Step2-1:潜在変数の変分事後分布の更新
				for(int n=0;n<m_numData;n++)
				{
					double gamma[MAXMODEL];

					// γの計算
					for(int m=0;m<=M;m++)
					{
						gamma[m]=CalcGamma(M,m,n);
					}

					// \bar{z}の計算
					for(int m=0;m<=M;m++)
					{
						double sum=0;
						for(int k=0;k<=M;k++)
						{
							sum+=exp(gamma[k]-gamma[m]);
						}
						params[m].z[n]=1.0/sum;
					}
				}

				//Step2-2:パラメータの変分事後分布の更新
				for(int m=0;m<=M;m++){
					PARAMETER &p = params[m];
					double invBeta[DIM][DIM] = {0};

					//初期化
					p.Nbar=0;
					for(int l=0;l<DIM;l++){
						p.xbar[l]=0;
						for(int ll=0;ll<DIM;ll++)
							p.C[l][ll]=0;
					}

					// \bar{x}_i, \bar{N}_i の計算
					for(int n=0;n<m_numData;n++){
						p.Nbar+=p.z[n];
						for(int l=0;l<DIM;l++)
							p.xbar[l]+=p.z[n]*m_data[l][n];
					}
					for(int l=0;l<DIM;l++)
						p.xbar[l]/=p.Nbar;

					// \bar{C}_i の計算
					for(int n=0;n<m_numData;n++){
						for(int l=0;l<DIM;l++)
							for(int ll=0;ll<DIM;ll++){
								p.C[l][ll] += p.z[n]*(m_data[l][n]-p.xbar[l])*(m_data[ll][n]-p.xbar[ll]);
							}
					}

					// φ, η, fの更新
					p.phi=m_phi0+p.Nbar;
					p.eta=m_eta0+p.Nbar;
					p.f=p.eta+1-DIM;

					//混合比
					p.alpha=p.Nbar/m_numData;

					// \bar{μ}, β, Σの更新
					InverseMatrix(p.beta,invBeta,DIM);	//+++
					for(int l=0;l<DIM;l++){
						p.mubar[l]=(p.Nbar*p.xbar[l]+m_xi0*m_nu0[m][l])/(p.Nbar+m_xi0);
						p.mu[l]=p.mubar[l];	//平均
						for(int ll=0;ll<DIM;ll++){
							p.beta[l][ll]=m_beta0[l][ll]+p.C[l][ll]+(p.Nbar*m_xi0)/(p.Nbar+m_xi0)*(p.xbar[l]-m_nu0[m][l])*(p.xbar[ll]-m_nu0[m][ll]);
							p.sigma[l][ll]=p.beta[l][ll]/((p.Nbar+m_xi0)*p.f);
							p.ivariance[l][ll]=p.eta*invBeta[l][ll];		//分散の逆行列(精度行列)(ウィシャート分布の期待値)
						}
					}
					InverseMatrix(p.ivariance,p.variance,DIM);			//分散(variance)行列
					for(int l=0;l<DIM;l++)
						for(int ll=0;ll<DIM;ll++)
							p.deviation[l][ll]=sqrt(p.variance[l][ll]);	//標準偏差行列
				}

				Fm = CalcFreeEnergy( m_models[M] , M );

				// 収束判定
				if( Fm - oldF < 0.00000001 )
				{
					break;
				}

				oldF = Fm;
			}

			// 自由エネルギーが最大の物を選択
			if( m_freeEnergies[M] < Fm )
			{
				m_freeEnergies[M] = Fm;
				SaveModelData( SaveDir , m_models[M],m_freeEnergies[M],M);	// 計算結果保存
			}
		}

		printf("計算完了...F = %lf\n" , m_freeEnergies[M] );
	}

	// 最適なモデルの選択
	double maxF = -DBL_MAX;
	for(int M=1 ; M<MAXMODEL ; M++ )
	{
		if( m_freeEnergies[M] > maxF )
		{
			maxF = m_freeEnergies[M];
			m_bestModelIdx = M;
		}
	}

}


// データをセット
int CVariationalBayes::SetData(double **data , int numData )
{
	//メモリ確保
	m_numData = numData;


	for(int i=0;i<DIM;i++)
	{
		if( m_data[i] ) free( m_data[i] );

		m_data[i] = new double[m_numData];
	}


	//データをコピー
	for(int cnt=0 ; cnt<m_numData ; cnt++)
	{
		for(int i=0;i<DIM;i++)
		{
			m_data[i][cnt] = data[cnt][i];
		}
	}

	return 0;
}


//初期パラメータの値決め
int CVariationalBayes::InitializeParameter(int model )
{
	m_phi0 = m_numData;
	m_xi0 = 1;
	m_eta0 = DIM + 2;

	double mu[DIM] = {0};	// 各次元の平均
	double var[DIM] = {0};	// 各次元の分散

	// 各次元の平均
	for(int n=0 ; n<m_numData ; n++ )
		for(int d=0 ; d<DIM ; d++ )	mu[d] += m_data[d][n];
	for(int d=0 ; d<DIM ; d++ )	mu[d] /= m_numData;

	// 各次元の分散
	for(int n=0 ; n<m_numData ; n++ )
		for(int d=0 ; d<DIM ; d++ )	var[d] += (m_data[d][n] - mu[d])*(m_data[d][n] - mu[d]);
	for(int d=0 ; d<DIM ; d++ )	var[d] /= m_numData;


	// ランダムでデータを選択
	for(int m=0;m<=model;m++)
	{
		int rnd = RandD() % m_numData;
		for(int d=0; d<DIM ; d++)
			m_nu0[m][d] = m_data[d][rnd];
	}


	// 分散
	for(int d1=0; d1<DIM ; d1++)
		for(int d0=0; d0<DIM ; d0++)
		{
			if( d0==d1 ) m_beta0[d1][d0] = var[d0];
			else m_beta0[d1][d0] = 0;
		}


	return 0;
}



//初期パラメータを代入する
void CVariationalBayes::SetInitializingParameter( PARAMETER* prm, int maxmodel, int model )
{
	prm->Nbar=(double)m_numData/(maxmodel+1);
	prm->phi=m_phi0 + prm->Nbar;
	prm->eta=m_eta0 + prm->Nbar;
	prm->f=m_eta0+prm->Nbar+1-DIM;

	for(int l=0;l<DIM;l++)
	{
		prm->mubar[l]=m_nu0[model][l];
		for(int ll=0;ll<DIM;ll++)
		{
			prm->beta[l][ll]=m_beta0[l][ll];
			prm->sigma[l][ll]=prm->beta[l][ll]/((prm->Nbar+m_xi0)*prm->f);
		}
	}
	return;
}

//γ[m][m_numData]を計算する
double CVariationalBayes::CalcGamma( int maxmodel, int model, int dataNo)
{
	double p1,p2,p3,p4,p5_1[DIM][DIM],p5_2[DIM][DIM],sum_p5[DIM][DIM],p5;
	double sunNbar,sumDigamma,iBeta[DIM][DIM],track;

	sunNbar = 0;
	sumDigamma = 0;

	for(int i=0;i<=maxmodel;i++)
		sunNbar+=m_models[maxmodel].m[i].Nbar;

	p1=digamma(m_phi0+m_models[maxmodel].m[model].Nbar);
	p2=-1*digamma((maxmodel+1)*m_phi0+sunNbar);

	for(int d=1;d<=DIM;d++)
		sumDigamma+=digamma((m_models[maxmodel].m[model].eta-d)/2);

	p3=0.5*sumDigamma;
	p4=(-0.5)*log(Determinant(m_models[maxmodel].m[model].beta,DIM));
	InverseMatrix(m_models[maxmodel].m[model].beta,iBeta,DIM);

	for(int i=0;i<DIM;i++)
	{
		for(int j=0;j<DIM;j++)
		{
			sum_p5[i][j]=0;
			p5_1[i][j]=m_models[maxmodel].m[model].eta*iBeta[i][j];
			p5_2[i][j]=(m_models[maxmodel].m[model].f/(m_models[maxmodel].m[model].f-2.0)*m_models[maxmodel].m[model].sigma[i][j])+(m_data[i][dataNo]-m_models[maxmodel].m[model].mubar[i])*(m_data[j][dataNo]-m_models[maxmodel].m[model].mubar[j]);
		}
	}

	for(int l=0;l<DIM;l++)
		for(int ll=0;ll<DIM;ll++)
			for(int k=0;k<DIM;k++)
				sum_p5[l][ll]+=p5_1[l][k]*p5_2[k][ll];

	track = TraceMatrix(sum_p5,DIM);
	p5 = (-0.5)*track;

	return p1+p2+p3+p4+p5;
}


//|A|:行列式
double CVariationalBayes::Determinant(double A[][DIM], int dim)
{
	double t,u,det;
	double S[DIM][DIM];

	for(int i=0;i<dim;i++)
		for(int j=0;j<dim;j++)
			S[i][j]=A[i][j];

	det = 1;
	for(int k=0;k<dim;k++){
		t = S[k][k];
		det *=t;
		for(int i=0;i<dim;i++) S[k][i] /=t;
		S[k][k] = 1/t;
		for(int j=0;j<dim;j++)
			if(j!=k){
				u=S[j][k];
				for(int i=0;i<dim;i++)
					if(i!=k) S[j][i] -= S[k][i] * u;
					else	 S[j][i] = -u/t;
			}
	}
	return det;
}

//逆行列:入力行列→A　出力用逆行列→S
void CVariationalBayes::InverseMatrix(double A[][DIM],double S[][DIM], int dim)
{
	double t,u,det;

	for(int i=0;i<dim;i++)
		for(int j=0;j<dim;j++)
			S[i][j]=A[i][j];

	det = 1;
	for(int k=0;k<dim;k++){
		t = S[k][k];
		det *=t;
		for(int i=0;i<dim;i++) S[k][i] /=t;
		S[k][k] = 1/t;
		for(int j=0;j<dim;j++)
			if(j!=k){
				u=S[j][k];
				for(int i=0;i<dim;i++)
					if(i!=k) S[j][i] -= S[k][i] * u;
					else	 S[j][i] = -u/t;
			}
	}
}

//行列のトレース::対角成分の和を計算
double CVariationalBayes::TraceMatrix(double S[][DIM], int dim)
{
	double sum=0;
	for(int i=0;i<dim;i++)
		sum+=S[i][i];
	return sum;
}

//自由エネルギーを計算する
double CVariationalBayes::CalcFreeEnergy( MODEL model, int m)
{
	double term[10] = {0};
	PARAMETER *p = model.m;
	double xi[MAXMODEL] = {0};
	double detB[MAXMODEL] = {0};
	double F = 0;

	// \xi_iの計算
	for(int i=0 ; i<=m ; i++ ) xi[i] = m_xi0 + p[i].Nbar;

	// |B_0|, |B_i|の計算
	double detB0 = Determinant( m_beta0 , DIM );
	for(int i=0 ; i<=m ; i++ ) detB[i] = Determinant( p[i].beta , DIM );

	// 第一項目
	for(int i=0 ; i<=m ; i++ )  term[0] += gamma( log( p[i].phi ) );

	// 第二項目
	term[1] = gamma( log( (m+1) * m_phi0) );

	// 第三項目
	term[2] = -gamma( log( (m+1)*m_phi0 + m_numData ) );

	// 第四項目
	term[3] = -gamma( log(m_phi0) ) * (m+1);

	// 第五項目
	term[4] = -( DIM * m_numData * log(M_PI)/log(2.0) )/2.0;

	// 第六項目
	for(int i=0 ; i<=m ; i++ ) term[5] += log( m_xi0 / xi[i] );
	term[5] *= DIM/2.0;

	// 第七項目
	term[6] = m_eta0 * log( detB0 ) / 2.0;
	for(int D=1 ; D<=DIM ; D++ ) term[6] -= gamma( log( (m_eta0+1-D) /2.0) );
	term[6] *= (m+1);

	// 第八項目
	for(int i=0 ; i<=m ; i++ )
	{
		term[7] += p[i].eta * log(detB[i]) /2.0;
		for(int D=1 ; D<=DIM ; D++ )
		{
			term[7] -= gamma( log( (p[i].eta+1-D)/2.0 ) );
		}
	}
	term[7] *= -1;

	// 総和
	for(int i=0 ; i<8 ; i++ )
	{
		F += term[i];
	}

	return F;
}

// 結果をファイルに保存する
void CVariationalBayes::SaveModelData( std::string saveDir , MODEL mod, double F , int modNum )
{
	double zmax,zresult;
	char strNo[32];
	char strBuff[256];
	sprintf( strNo , "%02d" , modNum+1 );
	FILE *fp = fopen( (saveDir + "result" + strNo + ".txt").c_str() , "w" );
	std::string filename = saveDir + "detail" + strNo + ".txt";


	// モデルのパラメータの保存
	SaveVariable( "■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■" , "設定値" ,  filename.c_str() , "w");
	SaveVariable( modNum+1 , "モデル数" , filename.c_str() );
	SaveVariable( m_numData , "データ数" , filename.c_str() );
	SaveVariable( F , "自由エネルギー" , filename.c_str() );

	SaveVariable( "■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■" , "\n\n初期値" , filename.c_str() );
	SaveVariable( m_phi0 , "φ" , filename.c_str() );
	SaveVariable( m_eta0 , "η" , filename.c_str() );
	SaveVariable( m_xi0 , "ξ" , filename.c_str() );
	SaveVariableMatrix( m_nu0 , modNum+1 , DIM , "ν" , filename.c_str() );

	for(int m=0;m<=modNum;m++)
	{
		sprintf( strBuff , "\n\nmodel = %d" , m );
		SaveVariable( "■■■■■■■■■■■■■■■■■■■■■■■■■■■■■■" , strBuff , filename.c_str() );

		sprintf( strBuff , "Nbar[%d]" , m );
		SaveVariable( mod.m[m].Nbar , strBuff , filename.c_str() );

		sprintf( strBuff , "phi[%d]" , m );
		SaveVariable( mod.m[m].phi , strBuff , filename.c_str() );

		sprintf( strBuff , "eta[%d]" , m );
		SaveVariable( mod.m[m].eta, strBuff , filename.c_str() );

		sprintf( strBuff , "f[%d]" , m );
		SaveVariable( mod.m[m].f , strBuff , filename.c_str() );

		sprintf( strBuff , "mubar[%d]" , m );
		SaveVariableArray( mod.m[m].mubar , DIM , strBuff , filename.c_str() );

		sprintf( strBuff , "xbar[%d]" , m );
		SaveVariableArray( mod.m[m].xbar , DIM , strBuff , filename.c_str() );

		sprintf( strBuff , "Beta[%d]" , m );
		SaveVariableMatrix( mod.m[m].beta , DIM , DIM , strBuff , filename.c_str() );

		sprintf( strBuff , "Sigma[%d]" , m );
		SaveVariableMatrix( mod.m[m].sigma , DIM , DIM , strBuff , filename.c_str() );

		sprintf( strBuff , "Cbar[%d]" , m );
		SaveVariableMatrix( mod.m[m].C , DIM , DIM , strBuff , filename.c_str() );

		sprintf( strBuff , "\nα[%d]" , m );
		SaveVariable( mod.m[m].alpha , strBuff , filename.c_str() );

		sprintf( strBuff , "平均[%d]" , m );
		SaveVariableArray( mod.m[m].mu , DIM , strBuff , filename.c_str() );

		sprintf( strBuff , "共分散[%d]" , m );
		SaveVariableMatrix( mod.m[m].variance , DIM , DIM , strBuff , filename.c_str() );

		sprintf( strBuff , "逆共分散[%d]" , m );
		SaveVariableMatrix( mod.m[m].ivariance , DIM , DIM , strBuff , filename.c_str() );
	}

	// クラスタリング結果を保存
	for(int n=0;n<m_numData;n++)
	{
		zresult=zmax=0;
		for(int m=0;m<=modNum;m++)
		{
			fprintf( fp ,"%.5lf	",mod.m[m].z[n]);
			if(zmax<mod.m[m].z[n])
			{
				zmax=mod.m[m].z[n];	
				zresult=m+1;	
			}
		}
		fprintf( fp ,"%.1lf\n",zresult);
	}
	fclose( fp );

}
