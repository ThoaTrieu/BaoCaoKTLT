#include <stdio.h>
#include <conio.h>
#include <windows.h>
#include <math.h>
#include <stdbool.h>
#include "giao_dien.h" 	
#include <graphics.h>
// cac bien toan cuc
const char* option[] =
{	
    "Tim cac KPL nghiem cua phuong trinh",
    "Thu hep khoang phan ly nghiem",
    "Ve do thi ham so f(x)",
    "Tim nghiem gan dung voi so lan lap n",
    "Tim nghiem gan dung voi sai so e",
    "Huong dan dieu khien MENU",
    "Gioi thieu ve phuong phap tiep tuyen",
    "Exit"
};

// ham in phuong trinh ra file
void print_poly_to_file (FILE *file){
	fprintf(file, "\n\nPhuong trinh: F(x) = %.*lfx^%d", decimal_digit, coeff[0], degree);
	for (int i = 1; i <= degree; i++){
    	if(coeff[i] >= 0)
        	fprintf(file, " +%.*lfx^%d", decimal_digit, coeff[i], degree - i);
        else 
        	fprintf(file, " %.*lfx^%d", decimal_digit, coeff[i], degree - i);
    }
    fprintf(file, " = 0\n");
}

// ham nhap phuong trinh
void input() {	
	int c;
    bool valid_input = false;
    FILE* file = fopen("input.txt", "w");
    
    if (file == NULL) {
    	TextColor(116);
        gotoxy(4,4); printf("%c Error. Khong the mo file de ghi.", 175);
        return;
    }
    
    while (!valid_input) {
    	khungngoai();
    	khungthongtin();
        TextColor(113); 
		gotoxy(4,4); printf("- Xet %c(x) la da thuc bac n co dang:", 159);
        gotoxy(4,6); printf(" %c(x) = an*x^n + a(n-1)*x^(n-1) + a(n-2)*x^(n-2) +...+ a0*x^0", 159);
        gotoxy(4,8); printf("Hay nhap da thuc %c(x) co dang nhu tren de bat dau tinh toan!", 159);
		gotoxy(4,10); printf("%c Nhap bac cua da thuc n = ", 175);
		
        if (scanf("%d", &degree) != 1 || degree < 0 || degree > N) {
        	TextColor(116);
        	gotoxy(4,12); printf("%c Error! Hay nhap mot so nguyen khong am nho hon 100.", 175);
        	gotoxy(4,13); printf("%c An Enter de nhap lai da thuc!", 175);
        	while (getchar() != '\n');  
            getchar();
            continue;
        }
        fprintf(file, "\nNguoi dung nhap vao: ");
        fprintf(file, "\n- Bac cua da thuc: %d.", degree);	
        
        TextColor(112); 
		gotoxy(4,12); printf("%c So chu so phan thap phan muon hien thi la: ", 175); 
		if (scanf("%d", &decimal_digit) != 1 || decimal_digit < 0 || decimal_digit > 20) {
        	TextColor(116);
        	gotoxy(4,13); printf("%c Error! Hay nhap mot so nguyen khong am nho hon 100.", 175);
        	gotoxy(4,14); printf("%c An Enter de nhap lai da thuc!", 175);
        	while (getchar() != '\n');
            getchar();
            continue;
        }
    	fprintf(file, "\n- So chu so phan thap phan muon hien thi la: %d.", decimal_digit);
    	
        TextColor(117);
		gotoxy(4,14); printf("%c Nhap cac he so cua da thuc: ", 175);
        fprintf(file, "\n\n- Cac he so cua da thuc:");
        for (int i = 0; i <= degree; i++) {
        	gotoxy(4,15 + i); printf("+) a%d = ", i);
            if (scanf("%lf", &coeff[i]) != 1) {
            	TextColor(116);
            	gotoxy(4,16 + i); printf("%c Error! He so khong hop le.", 175);
            	gotoxy(4,17 + i); printf("%c An Enter de nhap lai da thuc!", 175);
            	while (getchar() != '\n');  
            	getchar();
                valid_input = false;
                break;
            }
            fprintf(file, "\n +> a%d = %.*lf", i, decimal_digit, coeff[i]);
            valid_input = true;
        }
    }
    
    print_poly_to_file(file);
    fclose(file);
    
    TextColor(113);
    gotoxy(4,18 + degree); printf("%c Phuong trinh vua nhap la: ",175); 
    gotoxy(6,19 + degree); output(); 
    TextColor(116);
    gotoxy(4,21 + degree); printf("%c An Enter de bat dau su dung chuong trinh!",175);
    getch(); 
}

// ham in ra phuong trinh khi chay chuong trinh
void output()
{		
	printf("%c(x) = %.3lfx^%d", 159, coeff[0], degree);
    for (int i = 1; i <= degree; i++){
    	if(coeff[i] >= 0)
        	printf(" +%.3lfx^%d", coeff[i], degree - i);
        else 
        	printf(" %.3lfx^%d", coeff[i], degree - i);
    }
    printf(" = 0");
    TextColor(0);
}

// ham tinh gia tri f(x)
double f(double x)
{
    double temp = 0;
    for (int i = degree; i >= 0; i--) {
        temp += coeff[degree - i] * pow(x, i);
    }
    return temp;
}

// ham tinh gia tri dao ham cap 1 cua f(x)
double df(double x) {
    double h = 1e-7;
    return ((f(x + h) - f(x - h)) / (2 * h));
}

// ham tinh gia tri dao ham cap 2 cua f(x)
double ddf(double x) {
    double h = 1e-7;
    return ((df(x + h) - df(x - h)) / (2 * h));
}

// ham xac dinh dau cua mot gia tri x
double sign(double x)
{
	if (x>=0) return 1;
	else return -1;
}

// ham tra ve tri tuyet doi dao ham cap 1 cua f(x)
double f1(double x){
    return fabs(df(x));
}

// ham tra ve tri tuyet doi dao ham cap 2 cua f(x)
double f2(double x){
    return fabs(ddf(x));
}

// ham tinh sai so
double Compute_Error(double x, double x_old, double m1, double M2, int CT_SaiSo){
    if (CT_SaiSo == 1)
        return fabs(f(x))/ m1;
    else if (CT_SaiSo == 2)
        return M2 * pow((x - x_old), 2)/ (2 * m1);
}

// ham tim diem fourier
bool Fourier_Point(double a){
	if (sign(f(a)) == sign(ddf(a)))	
		return true;
    else return false;
}

// ham kiem tra dieu kien hoi tu
bool Check_Input(double a, double b){
    if (a == b) return false;
    else if (sign(f(a)) == sign(f(b))) return false;
    else if (a > b) return false;
    else if (sign(max(df, a, b)) != sign(min(df, a, b))) return false;
    else if (sign(max(ddf, a, b)) != sign(min(ddf, a, b))) return false;
    else return true;
}

// Ham tim ban kinh nghiem
double Find_Cauchy_Radius() {
    double maxRatio = 0;
    double an = coeff[0];
    for (int i = 1; i <= degree; i++) {
        if (coeff[i] != 0) {
            double ratio = fabs(coeff[i] / an);
            if (ratio > maxRatio) maxRatio = ratio;
        }
    }
    return 1 + maxRatio;
}

// Ham tim cac diem cuc tri 
int Find_Extrema_Points(double extrema[], double R, double h) {
    int count = 0;
    for (double x = -R; x < R; x += h) {
        double f1 = Derivative(x);
        double f2 = Derivative(x + h);
        if (f1 * f2 < 0 || fabs(f1) < EPSILON) {
            extrema[count++] = (x + x + h) / 2;
        }
    }
    return count;
}

//Ham tim cac khoang phan ly  nghiem 
int Find_Isolation_Intervals(double intervals[][2], double extrema[], int extremaCount, double R) {
    double points[extremaCount + 2];
    points[0] = -R;
    for (int i = 0; i < extremaCount; i++) points[i + 1] = extrema[i];
    points[extremaCount + 1] = R;

    int count = 0;
    for (int i = 0; i < extremaCount + 1; i++) {
        double a = points[i], b = points[i + 1];
        double fa = Evaluate(a), fb = Evaluate(b);
        if (fa * fb < 0) {
            intervals[count][0] = a;
            intervals[count++][1] = b;
        } else if (fabs(fa) < EPSILON) {
            intervals[count][0] = a - DELTA;
            intervals[count++][1] = a + DELTA;
        } else if (fabs(fb) < EPSILON) {
            intervals[count][0] = b - DELTA;
            intervals[count++][1] = b + DELTA;
        }
    }
    return count;
}

//Ham chon mpt khoang phan ly nghiem
void Select_Interval(double intervals[][2], int count, double *a, double *b) {
    *a = intervals[0][0];
    *b = intervals[0][1];
}

// Thu hep khoang bang phương phap chia doi
void Bisection_Interval(double *a, double *b) {
    while (fabs(*a - *b) > 0.5) {
        double mid = (*a + *b) / 2;
        if (Evaluate(*a) * Evaluate(mid) < 0) *b = mid;
        else *a = mid;
    }
}

// Kiem tra dieu kien hoi tu 
int Check_Input(double a, double b) {
    if (a == b) return 0;
    if (a > b) return 0;
    if (Evaluate(a) * Evaluate(b) > 0) return 0;
    return 1;
}

// Ham tim diem Fourier
int Fourier_Point(double x) {
    return (Evaluate(x) * SecondDerivative(x) > 0);
}

// Tinh X_n từ X_(n-1)
int Fourier_Point(double x) {
    return (Evaluate(x) * SecondDerivative(x) > 0);
}

// Ham tim min |f'(x)| trong [a,b]
double Min_Value(double a, double b) {
    double alpha = (b - a) / 10000;
    double minVal = fabs(Derivative(a));
    for (int i = 1; i <= 10000; i++) {
        double x = a + i * alpha;
        double fx = fabs(Derivative(x));
        if (fx < minVal) minVal = fx;
    }
    return minVal;
}

// Ham tim max |f''(x)| trong [a,b]
double Max_Value(double a, double b) {
    double alpha = (b - a) / 10000;
    double maxVal = fabs(SecondDerivative(a));
    for (int i = 1; i <= 10000; i++) {
        double x = a + i * alpha;
        double fx = fabs(SecondDerivative(x));
        if (fx > maxVal) maxVal = fx;
    }
    return maxVal;
}


// Tinh sai so 
void Compute_Error(double xn, double xn1, double m1, double M2, double *delta1, double *delta2) {
    *delta1 = fabs(Evaluate(xn)) / m1;
    *delta2 = (M2 / (2 * m1)) * (xn - xn1) * (xn - xn1);
}


// Vong lap chinh phuong phap Newton 
void Newton_Iteration(double x0, double m1, double M2) {
    double x1;
    int iteration = 0;
    double delta1, delta2;
    do {
        x1 = Newton_Step(x0);
        Compute_Error(x1, x0, m1, M2, &delta1, &delta2);
        printf("Iter %d: x = %.10lf, delta1 = %.10lf, delta2 = %.10lf\n", iteration, x1, delta1, delta2);
        x0 = x1;
        iteration++;
    } while ((delta1 >= EPSILON && delta2 >= EPSILON) && iteration < MAX_ITER);
    printf("Approximate root: %.10lf\n", x1);
}

int main() {
    printf("Nhap bac cua da thuc (<= %d): ", MAX_DEGREE);
    scanf("%d", &degree);
    printf("Nhap cac he so (tu bac %d den 0):\n", degree);
    for (int i = 0; i <= degree; i++) {
        scanf("%lf", &coeff[i]);
    }

    double R = Find_Cauchy_Radius();
    double extrema[100];
    int extremaCount = Find_Extrema_Points(extrema, R, 0.1);

    double intervals[100][2];
    int intervalCount = Find_Isolation_Intervals(intervals, extrema, extremaCount, R);

    double a, b;
    Select_Interval(intervals, intervalCount, &a, &b);
    Bisection_Interval(&a, &b);

    if (!Check_Input(a, b)) {
        printf("Khoang khong hop le.\n");
        return 1;
    }

    double x0 = Fourier_Point(a) ? a : b;
    double m1 = Min_Value(a, b);
    double M2 = Max_Value(a, b);

    Newton_Iteration(x0, m1, M2);
    return 0;
}

    