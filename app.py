import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import py3Dmol
from stmol import showmol
import time
import math

# --- دیکشنری نام‌ها (بدون تغییر) ---
smiles_to_name_map = {
    # ... (همان دیکشنری قبلی را اینجا کپی کنید) ...
    # C1-C7
    'C': 'Methane', 'CC': 'Ethane', 'CCC': 'Propane', 'CCCC': 'Butane', 'CC(C)C': 'Isobutane (2-Methylpropane)',
    'CCCCC': 'Pentane', 'CC(C)CC': 'Isopentane (2-Methylbutane)', 'CC(C)(C)C': 'Neopentane (2,2-Dimethylpropane)',
    'CCCCCC': 'Hexane', 'CC(C)CCC': 'Isohexane (2-Methylpentane)', 'CCC(C)CC': '3-Methylpentane',
    'CC(C)(C)CC': 'Neohexane (2,2-Dimethylbutane)', 'CC(C)C(C)C': '2,3-Dimethylbutane',
    'CCCCCCC': 'Heptane', 'CC(C)CCCC': '2-Methylhexane', 'CCC(C)CCC': '3-Methylhexane',
    'CC(C)C(C)CC': '2,3-Dimethylpentane', 'CC(C)CC(C)C': '2,4-Dimethylpentane',
    'CC(C)(C)CCC': '2,2-Dimethylpentane', 'CCC(C)(C)CC': '3,3-Dimethylpentane',
    'CC(CC)CCC': '3-Ethylpentane', 'C(C)(C)C(C)C': '2,2,3-Trimethylbutane', # یا CC(C)(C)C(C)C
    # C8 (partial)
    'CCCCCCCC': 'Octane', 'CC(C)CCCCC': '2-Methylheptane', 'CCC(C)CCCC': '3-Methylheptane', 'CCCC(C)CCC': '4-Methylheptane',
    'CC(CC)CCCC': '3-Ethylhexane','CC(C)(C)CCCC': '2,2-Dimethylhexane', # ... بقیه C8
    # C9 (partial)
    'CCCCCCCCC': 'Nonane', 'CC(C)CCCCCC': '2-Methylocatane', # ... بقیه C9
    # C10 (partial)
    'CCCCCCCCCC': 'Decane', 'CC(C)CCCCCCCC': '2-Methylnonane', # ... بقیه C10
}


# --- تابع گرفتن SMILES (بدون تغییر) ---
def get_alkane_isomer_smiles(n):
    # ... (همان تابع قبلی را اینجا کپی کنید) ...
    isomers_map = {
        1: [s for s in ['C'] if s in smiles_to_name_map], 2: [s for s in ['CC'] if s in smiles_to_name_map],
        3: [s for s in ['CCC'] if s in smiles_to_name_map], 4: [s for s in ['CCCC', 'CC(C)C'] if s in smiles_to_name_map],
        5: [s for s in ['CCCCC', 'CC(C)CC', 'CC(C)(C)C'] if s in smiles_to_name_map],
        6: [s for s in ['CCCCCC', 'CC(C)CCC', 'CCC(C)CC', 'CC(C)(C)CC', 'CC(C)C(C)C'] if s in smiles_to_name_map],
        7: [s for s in ['CCCCCCC', 'CC(C)CCCC', 'CCC(C)CCC', 'CC(C)C(C)CC', 'CC(C)CC(C)C','CC(C)(C)CCC', 'CCC(C)(C)CC', 'CC(CC)CCC', 'C(C)(C)C(C)C'] if s in smiles_to_name_map],
        8: [s for s in ['CCCCCCCC', 'CC(C)CCCCC', 'CCC(C)CCCC', 'CCCC(C)CCC','CC(CC)CCCC', 'CC(C)(C)CCCC', 'CCC(C)(C)CCC', 'CCCC(C)(C)CC','CC(C)C(C)CCC', 'CC(C)CC(C)CC', 'CC(C)CCC(C)C', 'CCC(C)C(C)CC','CC(CC)(C)CCC', 'CCC(C)(CC)CC', 'CC(C)(C)C(C)CC', 'CC(C)C(C)(C)CC','CC(C)(C)CC(C)C', 'CC(C)C(C)C(C)C', 'CC(C)(C)C(C)(C)C'] if s in smiles_to_name_map],
        9: [s for s in ['CCCCCCCCC', 'CC(C)CCCCCC', 'CCC(C)CCCCC', 'CCCC(C)CCCC', 'CCCCC(C)CCC','CC(CC)CCCCCC', 'CCC(CC)CCCCC', 'CCCC(CC)CCCC', 'CC(C)(C)CCCCCC','CCC(C)(C)CCCCC', 'CCCC(C)(C)CCCC', 'CCCCC(C)(C)CC', 'CC(C)C(C)CCCCC'] if s in smiles_to_name_map],
       10: [s for s in ['CCCCCCCCCC', 'CC(C)CCCCCCCC', 'CCC(C)CCCCCCC', 'CCCC(C)CCCCCC', 'CCCCC(C)CCCCC','CC(CC)CCCCCCCC', 'CCC(CC)CCCCCCC', 'CCCC(CC)CCCCCC', 'CCCCC(CC)CCCCC','CC(C)(C)CCCCCCCC', 'CCC(C)(C)CCCCCCC', 'CCCC(C)(C)CCCCCC', 'CCCCC(C)(C)CCCC','CCCCCC(C)(C)CC', 'CC(C)C(C)CCCCCCC'] if s in smiles_to_name_map]
    }
    return isomers_map.get(n, [])


# --- تابع ساخت *داده‌های* نمایش مولکول (نه خود نمایش) ---
def prepare_molecule_view_data(smiles_string, index):
    """Generates 3D coords and returns data needed for py3Dmol view."""
    try:
        mol = Chem.MolFromSmiles(smiles_string)
        if mol is None: return None # خطا را بعدا لاگ می‌کنیم
        mol = Chem.AddHs(mol)
        embed_result = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
        if embed_result == -1:
             try: AllChem.UFFOptimizeMolecule(mol); embed_result = AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
             except: pass
        # بهینه سازی اختیاری حتی اگر embed موفق بود
        if embed_result != -1:
             try: AllChem.UFFOptimizeMolecule(mol)
             except: pass # اگر بهینه‌سازی نشد مهم نیست

        # اگر embed ناموفق بود، None برگردان
        if embed_result == -1:
            st.warning(f"هشدار: تولید سه‌بعدی برای {smiles_string} ناموفق بود.")
            # اینجا می‌توانیم داده‌های نمایش دوبعدی را برگردانیم
            try:
                img = Draw.MolToImage(mol, size=(300,300))
                return {'type': '2d', 'image': img, 'smiles': smiles_string}
            except:
                return None # نمایش دوبعدی هم ممکن نبود
            # return None

        # اگر موفق بود، داده‌های سه‌بعدی را برگردان
        mol_block = Chem.MolToMolBlock(mol)
        return {'type': '3d', 'mol_block': mol_block, 'smiles': smiles_string}

    except Exception as e:
        st.error(f"خطای غیرمنتظره در آماده‌سازی {smiles_string}: {e}")
        return None

# --- تابع نمایش داده‌های مولکول ---
def display_molecule_data(view_data):
    """Displays the molecule based on prepared data."""
    if view_data is None:
        st.error("داده‌ای برای نمایش مولکول وجود ندارد.")
        return

    if view_data['type'] == '3d':
        try:
            view = py3Dmol.view(width=400, height=300)
            view.addModel(view_data['mol_block'], 'mol')
            view.setStyle({'stick': {'radius': 0.15}, 'sphere': {'scale': 0.25}})
            view.setBackgroundColor('#F5F5F5')
            view.zoomTo()
            showmol(view, height=400, width=400)
        except Exception as e:
            st.error(f"خطا در نمایش سه‌بعدی {view_data.get('smiles', '')}: {e}")
            # نمایش دوبعدی به عنوان جایگزین نهایی
            try:
                mol_2d = Chem.MolFromSmiles(view_data['smiles'])
                if mol_2d: img = Draw.MolToImage(mol_2d, size=(300,300)); st.image(img, caption=f"نمایش دوبعدی جایگزین: {view_data['smiles']}")
            except: pass

    elif view_data['type'] == '2d':
        try:
            st.image(view_data['image'], caption=f"نمایش دوبعدی: {view_data['smiles']}")
        except Exception as e:
             st.error(f"خطا در نمایش دوبعدی {view_data.get('smiles', '')}: {e}")
    else:
        st.error("نوع داده مولکول نامشخص است.")


# --- ساختار اصلی برنامه Streamlit ---

st.set_page_config(layout="wide", page_title="نمایشگر ایزومر آلکان")

st.title("🧪 نمایشگر ایزومرهای آلکان")
st.write("تعداد اتم‌های کربن (بین ۱ تا ۱۰) را وارد کنید تا ایزومرهای آن آلکان به صورت سه‌بعدی نمایش داده شوند.")
st.caption("توجه: لیست ایزومرها و نام‌ها برای تعداد کربن بالا کامل نیست.")

ITEMS_PER_PAGE = 5

if 'page_number' not in st.session_state: st.session_state.page_number = 0
if 'current_carbon_number' not in st.session_state: st.session_state.current_carbon_number = 5
if 'isomer_list' not in st.session_state: st.session_state.isomer_list = get_alkane_isomer_smiles(st.session_state.current_carbon_number)

with st.form("carbon_form"):
    carbon_number_input = st.number_input(
        label="تعداد کربن (n):", min_value=1, max_value=10,
        value=st.session_state.current_carbon_number, step=1,
        help="عددی بین ۱ تا ۱۰ وارد کنید."
    )
    submitted = st.form_submit_button(f"نمایش ایزومرهای C{carbon_number_input}H{2*carbon_number_input + 2}")
    if submitted:
        st.session_state.current_carbon_number = carbon_number_input
        st.session_state.isomer_list = get_alkane_isomer_smiles(carbon_number_input)
        st.session_state.page_number = 0

if st.session_state.isomer_list:
    all_isomer_smiles = st.session_state.isomer_list
    total_isomers = len(all_isomer_smiles)

    if total_isomers == 0:
         st.warning(f"برای n={st.session_state.current_carbon_number}، ایزومری در لیست یافت نشد.")
    else:
        st.success(f"تعداد {total_isomers} ایزومر برای C{st.session_state.current_carbon_number} یافت شد (نمایش {ITEMS_PER_PAGE} ایزومر در هر صفحه):")

        total_pages = math.ceil(total_isomers / ITEMS_PER_PAGE)
        if st.session_state.page_number >= total_pages and total_pages > 0: st.session_state.page_number = total_pages - 1
        elif st.session_state.page_number < 0: st.session_state.page_number = 0

        start_idx = st.session_state.page_number * ITEMS_PER_PAGE
        start_idx = max(0, start_idx)
        end_idx = min(start_idx + ITEMS_PER_PAGE, total_isomers)
        current_page_smiles = all_isomer_smiles[start_idx:end_idx]

        # --- مرحله جدید: آماده‌سازی داده‌های نمایش برای صفحه فعلی ---
        page_display_data = []
        with st.spinner("در حال آماده‌سازی مدل‌های سه‌بعدی..."):
            for i, smiles in enumerate(current_page_smiles):
                isomer_global_index = start_idx + i + 1
                view_data = prepare_molecule_view_data(smiles, isomer_global_index)
                if view_data: # فقط اگر آماده‌سازی موفق بود
                    isomer_name = smiles_to_name_map.get(smiles, "نام یافت نشد")
                    page_display_data.append({
                        'global_index': isomer_global_index,
                        'name': isomer_name,
                        'smiles': smiles,
                        'view_data': view_data
                    })
                else:
                    # اگر آماده‌سازی ناموفق بود، می‌توان یک placeholder اضافه کرد یا نادیده گرفت
                     page_display_data.append({
                        'global_index': isomer_global_index,
                        'name': smiles_to_name_map.get(smiles, "N/A"),
                        'smiles': smiles,
                        'view_data': None # علامت برای نمایش خطا
                    })


        # --- نمایش داده‌های آماده شده در ستون‌ها ---
        if page_display_data:
            num_columns = min(len(page_display_data), 3)
            if num_columns > 0:
                 cols = st.columns(num_columns)
                 # حالا بر اساس داده‌های آماده شده در ستون‌ها نمایش می‌دهیم
                 for i, display_item in enumerate(page_display_data):
                     col_index = i % num_columns
                     with cols[col_index]:
                         st.subheader(f"Isomer {display_item['global_index']} / {total_isomers}")
                         st.caption(f"**Name:** {display_item['name']}")
                         st.caption(f"`{display_item['smiles']}`")
                         with st.container():
                             # تابع جدید نمایش را صدا می‌زنیم
                             display_molecule_data(display_item['view_data'])


        st.markdown("---")

        # --- دکمه‌های ناوبری (بدون تغییر) ---
        col1, col2, col3 = st.columns([1, 2, 1])
        with col1:
            if st.button("⬅️ Previous", disabled=(st.session_state.page_number == 0)):
                st.session_state.page_number -= 1; st.rerun()
        with col2:
             if total_pages > 0: st.write(f"Page {st.session_state.page_number + 1} of {total_pages}")
             else: st.write("Page 0 of 0")
        with col3:
            if st.button("Next ➡️", disabled=(st.session_state.page_number >= total_pages - 1)):
                st.session_state.page_number += 1; st.rerun()

        st.markdown("---")
        st.info("💡 Use mouse to rotate, zoom, and pan the molecules.")

else:
    st.info("Please select the number of carbon atoms and click the button.")
