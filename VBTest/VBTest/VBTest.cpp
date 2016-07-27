/*******************************************************************

VBTest.cpp					developed by naka_t	2011.01.31

			�ϕ��x�C�YGMM�N���X�̎g�p��

Copyright (C) 2011  naka_t <naka_t@apple.ee.uec.ac.jp>
*******************************************************************/
#include "stdafx.h"
#include "VariationalBayes.h"
#include "utility.h"

int _tmain(int argc, _TCHAR* argv[])
{
	int dataNum = 0;
	int dim = 0;
	double **data = NULL;
	CVariationalBayes vb;

	// �f�[�^�̓ǂݍ���
	data = LoadMatrix<double>( dim , dataNum , "data.txt" );

	vb.SetData( data , dataNum );	// �f�[�^�����
	vb.Learn( "result" );			// �w�K�J�n

	printf("�œK������ %d\n" , vb.GetBestModelNum() );

	Free( data );
	return 0;
}

