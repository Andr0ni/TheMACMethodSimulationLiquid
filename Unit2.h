//---------------------------------------------------------------------------

#ifndef Unit2H
#define Unit2H
//---------------------------------------------------------------------------
#include <System.Classes.hpp>
#include <Vcl.Controls.hpp>
#include <Vcl.StdCtrls.hpp>
#include <Vcl.Forms.hpp>
#include <Vcl.ExtCtrls.hpp>
//---------------------------------------------------------------------------
class TForm2 : public TForm
{
__published:	// IDE-managed Components
    TButton *Button1;
    TImage *Image1;
    TLabeledEdit *LabeledEdit1;
    TLabeledEdit *LabeledEdit2;
    TLabeledEdit *LabeledEdit3;
    TLabeledEdit *LabeledEdit4;
    TLabeledEdit *LabeledEdit5;
    TColorBox *ColorBox1;
    TButton *Button2;
    TColorBox *ColorBox2;
    TCheckBox *CheckBox1;
    void __fastcall Button1Click(TObject *Sender);
    void __fastcall Button2Click(TObject *Sender);
    void __fastcall ColorBox1Change(TObject *Sender);
private:	// User declarations
public:		// User declarations
    __fastcall TForm2(TComponent* Owner);
};
//---------------------------------------------------------------------------
extern PACKAGE TForm2 *Form2;
//---------------------------------------------------------------------------
#endif
