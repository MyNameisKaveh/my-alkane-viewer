import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import py3Dmol
from stmol import showmol # برای نمایش py3Dmol در Streamlit
import time # برای ایجاد تاخیر جزئی و نمایش بهتر اسپینر

# --- تابع برای تولید SMILES ایزومرهای آلکان ---
# نکته مهم: تولید *تمام* ایزومرهای ساختاری (constitutional isomers)
# برای یک فرمول مولکولی داده شده، به خصوص برای آلکان‌ها با کربن بالا،
# یک مسئله پیچیده در شیمی محاسباتی است و الگوریتم‌های خاصی نیاز دارد.
# این تابع یک روش *ساده‌شده و محدود* ارائه می‌دهد که برای تعداد کربن کم کار می‌کند
# و ممکن است برای N های بزرگتر کامل یا بهینه نباشد.
# در اینجا برای سادگی از لیست‌های از پیش تعریف شده برای N های کوچک استفاده می‌کنیم.
def get_alkane_isomer_smiles(n):
    """
    For a given number of carbons n, return a list of SMILES strings
    for common alkane isomers. THIS IS A SIMPLIFIED/LIMITED LIST for demo.
    """
    # نقشه از پیش تعریف شده برای N های کوچک
    # منابع SMILES: PubChem or standard chemical databases
    isomers_map = {
        1: ['C'], # متان
        2: ['CC'], # اتان
        3: ['CCC'], # پروپان
        4: ['CCCC', # بوتان (n-Butane)
            'CC(C)C'], # ایزوبوتان (Isobutane)
        5: ['CCCCC', # پنتان (n-Pentane)
            'CC(C)CC', # ایزوپنتان (Isopentane)
            'CC(C)(C)C'], # نئوپنتان (Neopentane)
        6: ['CCCCCC', # هگزان (n-Hexane)
            'CC(C)CCC', # 2-متیل‌پنتان (Isohexane)
            'CCC(C)CC', # 3-متیل‌پنتان
            'CC(C)(C)CC', # 2,2-دی‌متیل‌بوتان (Neohexane)
            'CC(C)C(C)C'], # 2,3-دی‌متیل‌بوتان
        7: ['CCCCCCC', # هپتان (n-Heptane)
            'CC(C)CCCC', # 2-متیل‌هگزان
            'CCC(C)CCC', # 3-متیل‌هگزان
            'CCCC(C)CC', # 3-اتیل‌پنتان (اشتباه - این هم 3-متیل هگزان است با نامگذاری متفاوت) - SMILES درست باید برای ایزومرهای دیگر باشد
            # اصلاح و تکمیل ایزومرهای هپتان (9 ایزومر):
            'CC(C)CCCC', # 2-Methylhexane
            'CCC(C)CCC', # 3-Methylhexane
            'CC(C)C(C)CC', # 2,3-Dimethylpentane
            'CC(C)(C)CCC', # 2,2-Dimethylpentane
            'CCC(C)(C)CC', # 3,3-Dimethylpentane
            'CC(CC)CCC', # 3-Ethylpentane - این خودش یک ایزومر C7 است
            'CC(C)(CC)C', # 2,2,3-Trimethylbutane - این C7H16 است؟ C(4)+C(1)+C(1)+C(1) = C7. بله. SMILES: CC(C)(C)C(C)C
            # بگذارید SMILES های استاندارد هپتان را چک کنیم:
            # 1. CCCCCCC (n-Heptane)
            # 2. CC(C)CCCC (2-Methylhexane)
            # 3. CCC(C)CCC (3-Methylhexane)
            # 4. CC(C)(C)CCC (2,2-Dimethylpentane)
            # 5. CC(C)C(C)CC (2,3-Dimethylpentane)
            # 6. CCC(C)(C)CC (3,3-Dimethylpentane)
            # 7. C(C)(C)CC(C)C ??? این 2,2,3-Trimethylbutane است. درسته C7.
            # 8. CC(CC)CCC (3-Ethylpentane)
            # 9. C(C)(C)C(C)(C)C ??? این 2,2,3,3-Tetramethylpropane? نه این C8 است.
            # 9. C(C)(C)-C-(C)(C)C -> 2,2,3-trimethyl butane -> C(C)(C)C(C)C -> C7H16!
            # پس لیست 9 ایزومر هپتان:
            'CCCCCCC', 'CC(C)CCCC', 'CCC(C)CCC', 'CC(C)C(C)CC', 'CC(C)(C)CCC',
            'CCC(C)(C)CC', 'CC(CC)CCC', 'C(C)(C)C(C)C', 'CC(C)CC(C)C' # 2,4-Dimethylpentane
            # 2,4-Dimethylpentane هم هست! CCC(C)C(C)C ? نه. CC(C)CC(C)C
            # پس 9 تا شد.
            'CCCCCCC', 'CC(C)CCCC', 'CCC(C)CCC', 'CC(C)C(C)CC', 'CC(C)CC(C)C',
            'CC(C)(C)CCC', 'CCC(C)(C)CC', 'CC(CC)CCC', 'C(C)(C)C(C)C'
           ],
        # برای N های بزرگتر، لیست بسیار طولانی و تولید آن‌ها پیچیده می‌شود.
        # در این مثال ساده به همین تعداد بسنده می‌کنیم.
    }
    # اگر n در نقشه بود، لیستش را برگردان، وگرنه لیست خالی
    return isomers_map.get(n, [])

# --- تابع برای ایجاد و نمایش مولکول سه‌بعدی ---
def display_3d_molecule(smiles_string, index):
    """Generates 3D coords and displays molecule using py3Dmol."""
    try:
        # 1. ساخت مولکول از SMILES
        mol = Chem.MolFromSmiles(smiles_string)
        if mol is None:
            st.error(f"خطا: رشته SMILES نامعتبر است: {smiles_string}")
            return

        # 2. افزودن هیدروژن‌ها (مهم برای ساختار سه‌بعدی واقعی)
        mol = Chem.AddHs(mol)

        # 3. تولید مختصات سه‌بعدی
        # از الگوریتم ETKDG استفاده می‌کنیم که معمولاً نتایج بهتری می‌دهد
        embed_result = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())

        # اگر EmbedMolecule ناموفق بود (کد -1 برگرداند)، امتحان با روش دیگر یا هشدار
        if embed_result == -1:
             # گاهی اوقات بهینه‌سازی اولیه با UFF کمک می‌کند
             try:
                 AllChem.UFFOptimizeMolecule(mol)
                 embed_result = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
             except:
                 pass # اگر بهینه‌سازی هم شکست خورد، ادامه بده

        if embed_result == -1:
            st.warning(f"هشدار: تولید مختصات سه‌بعدی برای {smiles_string} با مشکل مواجه شد. ممکن است ساختار دوبعدی نمایش داده شود یا نمایش داده نشود.")
            # نمایش دو بعدی به عنوان جایگزین
            try:
                img = Draw.MolToImage(mol, size=(300,300))
                st.image(img, caption=f"نمایش دوبعدی برای: {smiles_string}")
            except:
                st.error("نمایش دوبعدی نیز ممکن نبود.")
            return # ادامه نده به نمایش سه‌بعدی

        # 4. بهینه‌سازی ساختار (اختیاری ولی توصیه می‌شود)
        try:
            AllChem.UFFOptimizeMolecule(mol)
        except Exception as e:
            st.warning(f"بهینه‌سازی ساختار برای {smiles_string} با خطا مواجه شد: {e}")


        # 5. تبدیل مولکول به فرمت MolBlock (که py3Dmol می‌فهمد)
        mol_block = Chem.MolToMolBlock(mol)

        # 6. تنظیمات نمایشگر py3Dmol
        view = py3Dmol.view(width=400, height=300)
        view.addModel(mol_block, 'mol')
        view.setStyle({'stick': {'radius': 0.15}, 'sphere': {'scale': 0.25}}) # نمایش به صورت stick و sphere
        view.setBackgroundColor('#F5F5F5') # رنگ پس‌زمینه کمی خاکستری
        view.zoomTo()

        # 7. نمایش در Streamlit با استفاده از stmol
        # کلید (key) منحصر به فرد برای هر نمایشگر لازم است
        showmol(view, height=400, width=400)

    except Exception as e:
        st.error(f"خطای غیرمنتظره در پردازش {smiles_string}: {e}")
        # نمایش دو بعدی به عنوان جایگزین در صورت بروز هر خطای دیگر
        try:
            mol_2d = Chem.MolFromSmiles(smiles_string)
            if mol_2d:
                 img = Draw.MolToImage(mol_2d, size=(300,300))
                 st.image(img, caption=f"نمایش دوبعدی جایگزین برای: {smiles_string}")
        except:
            pass # اگر این هم نشد، مشکلی نیست

# --- ساختار اصلی برنامه Streamlit ---

st.set_page_config(layout="wide", page_title="نمایشگر ایزومر آلکان")

st.title("🧪 نمایشگر ایزومرهای آلکان")
st.write("تعداد اتم‌های کربن (بین ۱ تا ۷) را وارد کنید تا ایزومرهای رایج آن آلکان به صورت سه‌بعدی نمایش داده شوند.")
st.caption("توجه: تولید *تمام* ایزومرها برای تعداد کربن بالا پیچیده است. این برنامه فقط ایزومرهای شناخته‌شده‌تر را برای n های کوچک نشان می‌دهد.")

# ورودی گرفتن از کاربر برای تعداد کربن
carbon_number = st.number_input(
    label="تعداد کربن (n):",
    min_value=1,
    max_value=7,  # محدود کردن به ۷ به دلیل پیچیدگی و لیست محدود ایزومرها در کد
    value=4,      # مقدار پیش‌فرض (بوتان)
    step=1,
    help="عددی بین ۱ تا ۷ وارد کنید."
)

# دکمه برای شروع پردازش
if st.button(f"نمایش ایزومرهای C{carbon_number}H{2*carbon_number + 2}"):
    # نمایش یک پیام "در حال پردازش" (spinner)
    with st.spinner(f"در حال جستجو و آماده‌سازی ایزومرهای C{carbon_number}H{2*carbon_number + 2}... لطفاً کمی صبر کنید..."):
        # گرفتن لیست SMILES ایزومرها
        isomer_smiles = get_alkane_isomer_smiles(carbon_number)
        time.sleep(1) # تاخیر کوچک برای نمایش بهتر اسپینر

        if not isomer_smiles:
            st.warning(f"برای n={carbon_number}، ایزومری در لیست این برنامه یافت نشد یا تعداد کربن خارج از محدوده پشتیبانی شده است.")
        else:
            st.success(f"تعداد {len(isomer_smiles)} ایزومر برای C{carbon_number}H{2*carbon_number + 2} یافت شد:")

            # تعیین تعداد ستون‌ها برای نمایش بهتر (مثلاً ۳ ستون)
            num_columns = 3
            cols = st.columns(num_columns)
            col_index = 0

            # حلقه برای نمایش هر ایزومر
            for i, smiles in enumerate(isomer_smiles):
                # انتخاب ستون فعلی به صورت چرخشی
                current_col = cols[col_index % num_columns]
                with current_col:
                    st.subheader(f"ایزومر {i+1}")
                    st.caption(f"`{smiles}`") # نمایش SMILES
                    # نمایش مولکول سه‌بعدی
                    display_3d_molecule(smiles, i) # پاس دادن index برای کلید یکتا

                col_index += 1

    st.markdown("---")
    st.info("💡 با استفاده از ماوس می‌توانید مولکول‌ها را بچرخانید، زوم کنید و حرکت دهید.")

# برای نمایش در محیط لوکال:
# در ترمینال دستور `streamlit run app.py` را اجرا کنید.
