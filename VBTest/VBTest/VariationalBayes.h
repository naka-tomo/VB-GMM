#pragma once
/*******************************************************************

VariationalBayes.h		developed by naka_t	2011.01.31

			変分ベイズGMM

　＊公開用にプログラムを整理					naka_t	2010.01.31

 Copyright (C) 2011  naka_t <naka_t@apple.ee.uec.ac.jp>
 *******************************************************************/
#include <string>

#define DIM			4		//データ次元数
#define MAXMODEL	10		//混合数の最大値
#define ITERATE		100		//1つの混合数に対して繰り返し計算する回数


class CVariationalBayes{
public:
	CVariationalBayes();
	~CVariationalBayes();

	void Learn( std::string SaveDir );													// 学習開始
	int SetData(double **data , int numData );											// データのセット
	double GetFreeEnergy(){ return m_freeEnergies[m_bestModelIdx]; }					// 最大の自由エネルギ－を返す
	double GetFreeEnergy(int m){ return m_freeEnergies[m-1]; }							// 特定のモデル数の自由エネルギーを返す
	int GetBestModelNum(){ return m_bestModelIdx+1; }									// 自由エネルギーが最大となった混合数を返す

protected:
	// モデルのパラメータ構造体
	typedef struct{
		double Nbar,phi,mubar[DIM],eta,beta[DIM][DIM],f,sigma[DIM][DIM];
		double *z,xbar[DIM],C[DIM][DIM],*gamma;
		double alpha,mu[DIM],variance[DIM][DIM],ivariance[DIM][DIM],deviation[DIM][DIM];
	}PARAMETER;

	// モデル構造体
	typedef struct{
		PARAMETER *m;
		double m_freeEnergies;
		int ModelNum;
	}MODEL;

	double *m_data[DIM];					// データ
	int m_numData;							// データ数
	int m_bestModelIdx;						// 最適なモデルのインデックス（モデル数+1）
	MODEL m_models[MAXMODEL];				// モデル

	double m_phi0,m_xi0,m_eta0,m_nu0[MAXMODEL][DIM],m_beta0[DIM][DIM];	// 初期パラメタ
	double m_freeEnergies[MAXMODEL];									// 自由エネルギー

	void FreeMemory();										// メモリ解放
	int AllocModel(MODEL *model , int modelnum );			// モデル確保
	int FreeModel( MODEL *model );							// モデル開放

	int InitializeParameter(int model );

	void SetInitializingParameter( PARAMETER* prm, int maxmodel, int model );
	double CalcGamma( int maxmodel, int model, int dataNo );
	double Determinant(double [][DIM],int dim);
	void InverseMatrix(double [][DIM],double [][DIM], int );
	double TraceMatrix(double[][DIM], int );
	double CalcFreeEnergy( MODEL mod, int modNum );
	void SaveModelData( std::string saveDir , MODEL mod, double F , int modNum );
};


