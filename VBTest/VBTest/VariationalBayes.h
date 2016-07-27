#pragma once
/*******************************************************************

VariationalBayes.h		developed by naka_t	2011.01.31

			�ϕ��x�C�YGMM

�@�����J�p�Ƀv���O�����𐮗�					naka_t	2010.01.31

 Copyright (C) 2011  naka_t <naka_t@apple.ee.uec.ac.jp>
 *******************************************************************/
#include <string>

#define DIM			4		//�f�[�^������
#define MAXMODEL	10		//�������̍ő�l
#define ITERATE		100		//1�̍������ɑ΂��ČJ��Ԃ��v�Z�����


class CVariationalBayes{
public:
	CVariationalBayes();
	~CVariationalBayes();

	void Learn( std::string SaveDir );													// �w�K�J�n
	int SetData(double **data , int numData );											// �f�[�^�̃Z�b�g
	double GetFreeEnergy(){ return m_freeEnergies[m_bestModelIdx]; }					// �ő�̎��R�G�l���M�|��Ԃ�
	double GetFreeEnergy(int m){ return m_freeEnergies[m-1]; }							// ����̃��f�����̎��R�G�l���M�[��Ԃ�
	int GetBestModelNum(){ return m_bestModelIdx+1; }									// ���R�G�l���M�[���ő�ƂȂ�����������Ԃ�

protected:
	// ���f���̃p�����[�^�\����
	typedef struct{
		double Nbar,phi,mubar[DIM],eta,beta[DIM][DIM],f,sigma[DIM][DIM];
		double *z,xbar[DIM],C[DIM][DIM],*gamma;
		double alpha,mu[DIM],variance[DIM][DIM],ivariance[DIM][DIM],deviation[DIM][DIM];
	}PARAMETER;

	// ���f���\����
	typedef struct{
		PARAMETER *m;
		double m_freeEnergies;
		int ModelNum;
	}MODEL;

	double *m_data[DIM];					// �f�[�^
	int m_numData;							// �f�[�^��
	int m_bestModelIdx;						// �œK�ȃ��f���̃C���f�b�N�X�i���f����+1�j
	MODEL m_models[MAXMODEL];				// ���f��

	double m_phi0,m_xi0,m_eta0,m_nu0[MAXMODEL][DIM],m_beta0[DIM][DIM];	// �����p�����^
	double m_freeEnergies[MAXMODEL];									// ���R�G�l���M�[

	void FreeMemory();										// ���������
	int AllocModel(MODEL *model , int modelnum );			// ���f���m��
	int FreeModel( MODEL *model );							// ���f���J��

	int InitializeParameter(int model );

	void SetInitializingParameter( PARAMETER* prm, int maxmodel, int model );
	double CalcGamma( int maxmodel, int model, int dataNo );
	double Determinant(double [][DIM],int dim);
	void InverseMatrix(double [][DIM],double [][DIM], int );
	double TraceMatrix(double[][DIM], int );
	double CalcFreeEnergy( MODEL mod, int modNum );
	void SaveModelData( std::string saveDir , MODEL mod, double F , int modNum );
};


