import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import py3Dmol
from stmol import showmol
import time
import math

# --- دیکشنری برای مپ کردن SMILES به نام آیوپاک ---
# منبع نام‌ها: PubChem / Wikipedia (نیاز به تکمیل و دقت دارد)
# نکته: پیدا کردن نام برای همه ۷۵ ایزومر دکان زمان‌بر است! فعلا چندتایی اضافه شده.
smiles_to_name_map = {
    # C1-C7
    'C': 'Methane',
    'CC': 'Ethane',
    'CCC': 'Propane',
    'CCCC': 'Butane',
    'CC(C)C': 'Isobutane (2-Methylpropane)',
    'CCCCC': 'Pentane',
    'CC(C)CC': 'Isopentane (2-Methylbutane)',
    'CC(C)(C)C': 'Neopentane (2,2-Dimethylpropane)',
    'CCCCCC': 'Hexane',
    'CC(C)CCC': 'Isohexane (2-Methylpentane)',
    'CCC(C)CC': '3-Methylpentane',
    'CC(C)(C)CC': 'Neohexane (2,2-Dimethylbutane)',
    'CC(C)C(C)C': '2,3-Dimethylbutane',
    'CCCCCCC': 'Heptane',
    'CC(C)CCCC': '2-Methylhexane',
    'CCC(C)CCC': '3-Methylhexane',
    'CC(C)C(C)CC': '2,3-Dimethylpentane',
    'CC(C)CC(C)C': '2,4-Dimethylpentane',
    'CC(C)(C)CCC': '2,2-Dimethylpentane',
    'CCC(C)(C)CC': '3,3-Dimethylpentane',
    'CC(CC)CCC': '3-Ethylpentane',
    'C(C)(C)C(C)C': '2,2,3-Trimethylbutane', # یا CC(C)(C)C(C)C
    # C8 (18 isomers - partial names for demo)
    'CCCCCCCC': 'Octane',
    'CC(C)CCCCC': '2-Methylheptane',
    'CCC(C)CCCC': '3-Methylheptane',
    'CCCC(C)CCC': '4-Methylheptane',
    'CC(CC)CCCC': '3-Ethylhexane',
    'CC(C)(C)CCCC': '2,2-Dimethylhexane',
    'CCC(C)(C)CCC': '3,3-Dimethylhexane',
    'CCCC(C)(C)CC': 'Name Needed', # 3,4-Dimethylhexane
    'CC(C)C(C)CCC': '2,3-Dimethylhexane',
    'CC(C)CC(C)CC': '2,4-Dimethylhexane',
    'CC(C)CCC(C)C': '2,5-Dimethylhexane',
    'CCC(C)C(C)CC': '3,4-Dimethylhexane',
    'CC(CC)(C)CCC': '3-Ethyl-2-methylpentane', # Check IUPAC vs SMILES
    'CCC(C)(CC)CC': '3-Ethyl-3-methylpentane',
    'CC(C)(C)C(C)CC': '2,2,3-Trimethylpentane',
    'CC(C)C(C)(C)CC': '2,3,3-Trimethylpentane',
    'CC(C)(C)CC(C)C': '2,2,4-Trimethylpentane',
    'CC(C)C(C)C(C)C': '2,3,4-Trimethylpentane',
    'CC(C)(C)C(C)(C)C': '2,2,3,3-Tetramethylbutane',
    # C9 (35 isomers - only first few named)
    'CCCCCCCCC': 'Nonane',
    'CC(C)CCCCCC': '2-Methylocatane',
    'CCC(C)CCCCC': '3-Methylocatane',
    'CCCC(C)CCCC': '4-Methylocatane',
    'CCCCC(C)CCC': '5-Methylocatane', # same as 4-Methyloctane due to symmetry
    'CC(CC)CCCCCC': '3-Ethylheptane',
    'CCC(CC)CCCCC': '4-Ethylheptane',
    'CCCC(CC)CCCC': 'Name Needed', # 4-Propylhexane? No. Ethyl+Methyl?
    'CC(C)(C)CCCCCC': '2,2-Dimethylheptane',
    'CCC(C)(C)CCCCC': '3,3-Dimethylheptane',
    'CCCC(C)(C)CCCC': '4,4-Dimethylheptane',
    'CCCCC(C)(C)CC': 'Name Needed',
    'CC(C)C(C)CCCCC': '2,3-Dimethylheptane',
    # C10 (75 isomers - only first few named)
    'CCCCCCCCCC': 'Decane',
    'CC(C)CCCCCCCC': '2-Methylnonane',
    'CCC(C)CCCCCCC': '3-Methylnonane',
    'CCCC(C)CCCCCC': '4-Methylnonane',
    'CCCCC(C)CCCCC': '5-Methylnonane',
    'CC(CC)CCCCCCCC': '3-Ethyloctane',
    'CCC(CC)CCCCCCC': '4-Ethyloctane',
    'CCCC(CC)CCCCCC': 'Name Needed',
    'CCCCC(CC)CCCCC': 'Name Needed',
    'CC(C)(C)CCCCCCCC': '2,2-Dimethylocatane',
    'CCC(C)(C)CCCCCCC': '3,3-Dimethylocatane',
    'CCCC(C)(C)CCCCCC': '4,4-Dimethylocatane',
    'CCCCC(C)(C)CCCC': '5,5-Dimethylocatane', # Same as 4,4? No.
    'CCCCCC(C)(C)CC': 'Name Needed',
    'CC(C)C(C)CCCCCCC': '2,3-Dimethylocatane',
}

# --- تابع برای گرفتن SMILES (از کلیدهای مپ بالا استفاده می‌کنیم تا مطمئن باشیم نام وجود دارد) ---
def get_alkane_isomer_smiles(n):
    """
    Returns a list of SMILES strings for alkane isomers based on available names.
    """
    # فیلتر کردن دیکشنری نام‌ها بر اساس تعداد کربن (تقریبی و نه دقیق)
    # راه بهتر: از خود SMILES تعداد کربن را بشماریم
    smiles_list = []
    for smiles in smiles_to_name_map.keys():
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol and mol.GetNumAtoms(onlyExplicit=False) - mol.GetNumHeavyAtoms() == n: # GetNumHeavyAtoms gives non-H count (approx C count for alkanes)
                 # Check exact formula CnH(2n+2)
                 formula = Chem.rdMolDescriptors.CalcMolFormula(mol)
                 expected_formula = f"C{n}H{2*n+2}"
                 # Sometimes formula might have implicit H like C10H22. Check just carbon count.
                 atom_counts = {}
                 for atom in mol.GetAtoms():
                     symbol = atom.GetSymbol()
                     atom_counts[symbol] = atom_counts.get(symbol, 0) + 1
                 if atom_counts.get('C', 0) == n:
                      smiles_list.append(smiles)
        except:
            continue # Skip if SMILES is invalid

    # روش ساده‌تر: بر اساس پیشوند فرمول در PubChem (اگر از آنجا کپی شده باشد)
    # یا فقط بر اساس تعریف دستی لیست‌ها مثل قبل (که کنترل بیشتری می‌دهد)
    # بگذارید به روش قبلی (لیست دستی) برگردیم چون دقیق‌تر است

    isomers_map = {
        1: [s for s in ['C'] if s in smiles_to_name_map],
        2: [s for s in ['CC'] if s in smiles_to_name_map],
        3: [s for s in ['CCC'] if s in smiles_to_name_map],
        4: [s for s in ['CCCC', 'CC(C)C'] if s in smiles_to_name_map],
        5: [s for s in ['CCCCC', 'CC(C)CC', 'CC(C)(C)C'] if s in smiles_to_name_map],
        6: [s for s in ['CCCCCC', 'CC(C)CCC', 'CCC(C)CC', 'CC(C)(C)CC', 'CC(C)C(C)C'] if s in smiles_to_name_map],
        7: [s for s in ['CCCCCCC', 'CC(C)CCCC', 'CCC(C)CCC', 'CC(C)C(C)CC', 'CC(C)CC(C)C',
                        'CC(C)(C)CCC', 'CCC(C)(C)CC', 'CC(CC)CCC', 'C(C)(C)C(C)C'] if s in smiles_to_name_map],
        8: [s for s in ['CCCCCCCC', 'CC(C)CCCCC', 'CCC(C)CCCC', 'CCCC(C)CCC',
                        'CC(CC)CCCC', 'CC(C)(C)CCCC', 'CCC(C)(C)CCC', 'CCCC(C)(C)CC', # Missing SMILES for one name?
                        'CC(C)C(C)CCC', 'CC(C)CC(C)CC', 'CC(C)CCC(C)C', 'CCC(C)C(C)CC',
                        'CC(CC)(C)CCC', 'CCC(C)(CC)CC', 'CC(C)(C)C(C)CC', 'CC(C)C(C)(C)CC',
                        'CC(C)(C)CC(C)C', 'CC(C)C(C)C(C)C', 'CC(C)(C)C(C)(C)C'] if s in smiles_to_name_map],
        9: [s for s in ['CCCCCCCCC', 'CC(C)CCCCCC', 'CCC(C)CCCCC', 'CCCC(C)CCCC', 'CCCCC(C)CCC',
                        'CC(CC)CCCCCC', 'CCC(CC)CCCCC', 'CCCC(CC)CCCC', 'CC(C)(C)CCCCCC',
                        'CCC(C)(C)CCCCC', 'CCCC(C)(C)CCCC', 'CCCCC(C)(C)CC', 'CC(C)C(C)CCCCC'] if s in smiles_to_name_map],
       10: [s for s in ['CCCCCCCCCC', 'CC(C)CCCCCCCC', 'CCC(C)CCCCCCC', 'CCCC(C)CCCCCC', 'CCCCC(C)CCCCC',
                        'CC(CC)CCCCCCCC', 'CCC(CC)CCCCCCC', 'CCCC(CC)CCCCCC', 'CCCCC(CC)CCCCC',
                        'CC(C)(C)CCCCCCCC', 'CCC(C)(C)CCCCCCC', 'CCCC(C)(C)CCCCCC', 'CCCCC(C)(C)CCCC',
                        'CCCCCC(C)(C)CC', 'CC(C)C(C)CCCCCCC'] if s in smiles_to_name_map]
    }

    return isomers_map.get(n, [])


# --- تابع نمایش مولکول (بدون تغییر) ---
def display_3d_molecule(smiles_string, index):
    """Generates 3D coords and displays molecule using py3Dmol."""
    try:
        mol = Chem.MolFromSmiles(smiles_string)
        if mol is None: st.error(f"خطا: SMILES نامعتبر: {smiles_string}"); return
        mol = Chem.AddHs(mol)
        embed_result = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        if embed_result == -1:
             try: AllChem.UFFOptimizeMolecule(mol); embed_result = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
             except: pass
        if embed_result == -1:
            st.warning(f"هشدار: تولید مختصات سه‌بعدی برای {smiles_string} ناموفق بود.")
            try: img = Draw.MolToImage(mol, size=(300,300)); st.image(img, caption=f"نمایش دوبعدی: {smiles_string}")
            except: st.error("نمایش دوبعدی نیز ممکن نبود."); return
            return
        try: AllChem.UFFOptimizeMolecule(mol)
        except Exception as e: st.warning(f"بهینه‌سازی ساختار برای {smiles_string} ناموفق: {e}")

        mol_block = Chem.MolToMolBlock(mol)
        view = py3Dmol.view(width=400, height=300)
        view.addModel(mol_block, 'mol')
        view.setStyle({'stick': {'radius': 0.15}, 'sphere': {'scale': 0.25}})
        view.setBackgroundColor('#F5F5F5')
        view.zoomTo()
        # Pass a unique key based on index and smiles for stability within loops/reruns
        # Note: Even though showmol doesn't take 'key', wrapping it might help Streamlit
        # Let's try without key first, as per previous fix. If instability occurs, consider st.container with key.
        showmol(view, height=400, width=400)

    except Exception as e:
        st.error(f"خطای غیرمنتظره در پردازش {smiles_string}: {e}")
        try:
            mol_2d = Chem.MolFromSmiles(smiles_string)
            if mol_2d: img = Draw.MolToImage(mol_2d, size=(300,300)); st.image(img, caption=f"نمایش دوبعدی جایگزین: {smiles_string}")
        except: pass

# --- ساختار اصلی برنامه Streamlit ---

st.set_page_config(layout="wide", page_title="نمایشگر ایزومر آلکان")

st.title("🧪 نمایشگر ایزومرهای آلکان")
st.write("تعداد اتم‌های کربن (بین ۱ تا ۱۰) را وارد کنید تا ایزومرهای آن آلکان به صورت سه‌بعدی نمایش داده شوند.")
st.caption("توجه: لیست ایزومرها و نام‌ها برای تعداد کربن بالا کامل نیست.")

# --- پارامترهای صفحه‌بندی ---
ITEMS_PER_PAGE = 5

# --- مقداردهی اولیه Session State ---
if 'page_number' not in st.session_state:
    st.session_state.page_number = 0
if 'current_carbon_number' not in st.session_state:
    st.session_state.current_carbon_number = 5 # مقدار پیش‌فرض اولیه
if 'isomer_list' not in st.session_state:
    # لیست اولیه را بر اساس کربن پیش‌فرض پر می‌کنیم
    st.session_state.isomer_list = get_alkane_isomer_smiles(st.session_state.current_carbon_number)

# --- ورودی گرفتن از کاربر با فرم ---
with st.form("carbon_form"):
    # مقدار پیش‌فرض number_input را از session_state می‌خوانیم
    carbon_number_input = st.number_input(
        label="تعداد کربن (n):",
        min_value=1,
        max_value=10,
        value=st.session_state.current_carbon_number, # خواندن از state
        step=1,
        help="عددی بین ۱ تا ۱۰ وارد کنید."
    )
    submitted = st.form_submit_button(f"نمایش ایزومرهای C{carbon_number_input}H{2*carbon_number_input + 2}")

    # وقتی فرم سابمیت شد، state را آپدیت می‌کنیم
    if submitted:
        st.session_state.current_carbon_number = carbon_number_input
        st.session_state.isomer_list = get_alkane_isomer_smiles(carbon_number_input)
        st.session_state.page_number = 0 # ریست کردن صفحه به اول
        # نیازی به st.rerun() نیست چون Streamlit بعد از submit فرم خودش اجرا می‌شود

# --- نمایش ایزومرها (خارج از فرم) ---
# فقط اگر لیستی برای نمایش در state وجود داشته باشد
if st.session_state.isomer_list:
    all_isomer_smiles = st.session_state.isomer_list
    total_isomers = len(all_isomer_smiles)

    if total_isomers == 0:
         st.warning(f"برای n={st.session_state.current_carbon_number}، ایزومری در لیست یافت نشد.")
    else:
        st.success(f"تعداد {total_isomers} ایزومر برای C{st.session_state.current_carbon_number} یافت شد (نمایش {ITEMS_PER_PAGE} ایزومر در هر صفحه):")

        # محاسبه تعداد کل صفحات
        total_pages = math.ceil(total_isomers / ITEMS_PER_PAGE)

        # اطمینان از اینکه شماره صفحه معتبر است
        if st.session_state.page_number >= total_pages and total_pages > 0:
            st.session_state.page_number = total_pages - 1
        elif st.session_state.page_number < 0:
            st.session_state.page_number = 0


        # محاسبه اندیس شروع و پایان برای صفحه فعلی
        start_idx = st.session_state.page_number * ITEMS_PER_PAGE
        # اطمینان از اینکه اندیس شروع منفی یا خیلی بزرگ نشود
        start_idx = max(0, start_idx)
        end_idx = min(start_idx + ITEMS_PER_PAGE, total_isomers) # جلوگیری از بیرون زدن از لیست

        # گرفتن لیست ایزومرهای صفحه فعلی
        current_page_isomers = all_isomer_smiles[start_idx:end_idx]

        # نمایش ایزومرهای صفحه فعلی در ستون‌ها
        if current_page_isomers: # فقط اگر ایزومری برای نمایش در این صفحه هست
            num_columns = min(len(current_page_isomers), 3) # ستون‌ها بر اساس تعداد واقعی آیتم‌ها
            if num_columns > 0: # اگر حداقل یک ایزومر هست
                 cols = st.columns(num_columns)
                 for i, smiles in enumerate(current_page_isomers):
                     col_index = i % num_columns
                     with cols[col_index]:
                         # نمایش اندیس کلی ایزومر
                         isomer_global_index = start_idx + i + 1
                         isomer_name = smiles_to_name_map.get(smiles, "نام یافت نشد") # گرفتن نام

                         st.subheader(f"Isomer {isomer_global_index} / {total_isomers}")
                         st.caption(f"**Name:** {isomer_name}") # نمایش نام
                         st.caption(f"`{smiles}`") # نمایش SMILES
                         # استفاده از container با key برای پایداری بهتر ویجت‌ها در rerun
                         with st.container():
                              display_3d_molecule(smiles, isomer_global_index)

        st.markdown("---") # خط جداکننده

        # --- دکمه‌های ناوبری صفحه‌بندی ---
        col1, col2, col3 = st.columns([1, 2, 1])

        with col1:
            # دکمه "قبلی"
            if st.button("⬅️ Previous", disabled=(st.session_state.page_number == 0)):
                st.session_state.page_number -= 1
                st.rerun() # صفحه را دوباره بارگذاری کن

        with col2:
             if total_pages > 0:
                  st.write(f"Page {st.session_state.page_number + 1} of {total_pages}")
             else:
                  st.write("Page 0 of 0")


        with col3:
            # دکمه "بعدی"
            if st.button("Next ➡️", disabled=(st.session_state.page_number >= total_pages - 1)):
                st.session_state.page_number += 1
                st.rerun() # صفحه را دوباره بارگذاری کن

        st.markdown("---")
        st.info("💡 Use mouse to rotate, zoom, and pan the molecules.")

else:
    # اگر state خالی بود (مثلاً در اولین اجرا قبل از سابمیت)
    st.info("Please select the number of carbon atoms and click the button.")
