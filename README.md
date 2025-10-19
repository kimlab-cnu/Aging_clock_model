# Aging-Clock-Model
- 노화연구는 단백질 발현량 변화 쪽으로 많이 연구되어 왔지만, 노화는 단순히 단백질 발현량 변화만으로 설명되지 않습니다.
- 저희는 특히 아직 충분히 연구되지 않은 영역인 단백질 구조적 변화와 상호작용 네트워크의 dynamics에 주목했습니다.
- 이번 연구에서는 혈장 단백질의 age dependent한 structural dynamics를 포괄적으로 분석하기 위해 실험적 프로테오믹스와 computational 구조생물학적 접근을 통합하여, **혈장 단백질이 건강한 삶의 주기에 따라 어떻게 구조적으로 변화하는지**를 전반적으로 분석했습니다. 


---
## 1. Data preprocessing
- #### Data_preprocessing_for_PAR.R
:Normalize에 사용할 internal standard(IS)를 선택하기 위해서, IS 중 peptide-03의 1,2,3번째 precursor의 6개 이온을 이용하여 각 이온의 PAR (Peak Area Ratio) 값을 계산함.

(251019 피규어1B 코드 추가 완료)


- #### PAR_distribution_histogram.R
:위에서 구한 PAR1~6에 log10을 씌우고, CV값도 같이 나타낸 히스토그램을 그림.

(251019 피규어1A,C 코드 추가 완료)


- #### PAR_input_preprocessing.R
- #### Aging_log_norm_data_preprocessing.ipynb
- #### PAR_data_preprocessing_Jhyh.R
:이후 분석에서 필요한 데이터를 전처리.


## 2. Statistical Analysis

- #### Data_preprocessing_Jhyh+cosine similatiry.R
:전처리한 데이터 파일을 사용하여 코사인 유사도 분석을 진행.

(피규어2A,B)


- #### limma_test_Jh_수정_추가_fit3_최종.R
:limma를 사용하여 패턴별 FC, p-value 조건을 충족하는 이온만 추출, 패턴별 평균 PAR 변화 그래프를 그림.

(피규어3A,B)

- #### overlap_pattern_fit3_최종.R
:코사인 유사도 분석 결과와 limma 결과를 합쳐, 중복되는 이온만 추출.


- #### overlap_boxplot_fit3_최종.R
:패턴별 각 이온의 박스플랏을 그림.

(피규어3C)

## 3. Structural Analysis

- #### pymol 상호작용 거리 계산 코드.txt
:AlphaFold3 multimer 결과를 pymol로 열어서 확인할 때, 패턴별 단백질 복합체(1:1) 내 펩타이드 및 단백질 인터페이스의 최단거리 확인함.


- #### pymol 상호작용 거리 계산 코드.txt
:AlphaFold3 multimer 결과를 pymol로 열어서 확인할 때, mean pLDDT를 구함

(피규어4B,C,D)

---
- ##### 참고
cytoscape: Full STRING Network / confidence ≥ 0.7 / Homo sapiens

AF3 multimer 웹서버: PTM은 N-linked glycosylation site만 NAG로 추가

