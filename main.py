import streamlit as st
import pandas as pd
# from streamlit_molstar import st_molstar, st_molstar_rcsb, st_molstar_remote
# 如需使用口袋预测相关函数
from streamlit_molstar.pocket import (
    select_pocket_from_local_protein,
    # 如果你的项目需要也可以 import select_pocket_from_upload_protein
)
# docking 模块
from streamlit_molstar.docking import st_molstar_docking

import os
import json
import subprocess
import tempfile  # 用于创建临时文件
import re
# import tqdm

# 如果没有在 session_state 中记录 page,就初始化一个默认值
if 'page' not in st.session_state:
    st.session_state['page'] = '主页'

# 在侧边栏使用按钮来切换页面
st.sidebar.title("Navigation")
if st.sidebar.button("主页"):
    st.session_state['page'] = "主页"

# 新增“准备配体”按钮（插在“主页”和“口袋预测”之间）
if st.sidebar.button("准备配体"):
    st.session_state['page'] = "准备配体"

if st.sidebar.button("口袋预测"):
    st.session_state['page'] = "口袋预测"
if st.sidebar.button("分子对接"):
    st.session_state['page'] = "分子对接"
if st.sidebar.button("批量口袋预测与对接"):
    st.session_state['page'] = "批量口袋预测与对接"
    # 新增“预测亲和力”按钮
if st.sidebar.button("预测亲和力"):
    st.session_state['page'] = "预测亲和力"
# 新增“作者信息”按钮
if st.sidebar.button("作者信息"):
    st.session_state['page'] = "作者信息"

# 获取当前页面
page = st.session_state['page']

# ------------------------------------------------------------------------------
# 主页
# ------------------------------------------------------------------------------
if page == "主页":
    # 使用 HTML 和 Markdown 居中标题
    st.markdown(
        "<h1 style='text-align: center;'>⚡️欢迎使用⚡️</h1>",
        unsafe_allow_html=True,
    )
    st.markdown("<br>", unsafe_allow_html=True)

    # 显示字符画
    try:
        with open("./lib/logo.txt", "r", encoding="utf-8") as file:
            ascii_art = file.read()
            styled_ascii_art = ascii_art.replace(" ", "&nbsp;").replace("\n", "<br>")
            html_code = f"""
            <div style='text-align: center; font-family: monospace; font-size: 14px; line-height: 1;'>
                {styled_ascii_art}
            </div>
            """
            st.markdown(html_code, unsafe_allow_html=True)

    except FileNotFoundError:
        st.error("logo.txt 文件未找到,请确保它与脚本位于同一目录下")
    except UnicodeDecodeError:
        st.error("无法解码 logo.txt 文件,请确认文件编码格式是否为 UTF-8")

    # 在字符画和图片之间插入若干空行
    st.markdown("<br><br><br><br><br>", unsafe_allow_html=True)

    # 显示 logo.png
    if os.path.exists("./lib/logo.png"):
        st.image("./lib/logo.png", use_container_width=True)
    else:
        st.error("logo.png 文件未找到,请确保它位于 lib 目录下")

# ------------------------------------------------------------------------------
# 准备配体
# ------------------------------------------------------------------------------
elif page == "准备配体":
    st.title("准备配体")

    import os
    from rdkit import Chem
    from rdkit.Chem import AllChem, Draw
    import py3Dmol
    from stmol import showmol
    from streamlit_ketcher import st_ketcher

    # 1. 允许用户上传一个 SDF 文件
    st.info("请上传一个 SDF 文件,或在画布中绘制分子结构/粘贴SMILES")
    st.markdown("**上传**分子文件 (SDF 格式): ")
    sdf_file = st.file_uploader("", type=["sdf"])
    # 2. 允许用户使用 Ketcher 绘制或输入 SMILES
    st.markdown("**或者** 在下方绘制分子结构/粘贴SMILES: ")
    smiles_input = st_ketcher()

    def process_and_show_mol(
        mol: Chem.Mol, 
        uploaded_sdf_name: str = None, 
        user_defined_filename: str = None
    ):
        """
        对分子进行加氢、3D 嵌入、MMFF 优化并展示 2D/3D 结构；
        根据不同来源决定最终保存的文件名：
        - 如果有 uploaded_sdf_name, 则用 "原文件名去除.sdf + '_prepared.sdf'"
        - 如果没有 uploaded_sdf_name, 但用户给了自定义文件名, 则用 "用户自定义文件名 + '.sdf'"
        """
        if not mol:
            return

        # 2D 可视化
        st.subheader("2D 分子结构")
        st.image(Draw.MolToImage(mol, size=(300, 300)), use_container_width=False)

        # 生成 3D 构象并能量优化
        mol_3d = Chem.AddHs(mol)
        AllChem.EmbedMolecule(mol_3d, AllChem.ETKDG())
        AllChem.MMFFOptimizeMolecule(mol_3d)

        # 3D 可视化
        st.subheader("3D 分子结构")
        mol_block = Chem.MolToMolBlock(mol_3d)
        xyzview = py3Dmol.view(width=500, height=400)
        xyzview.addModel(mol_block, "mol")
        xyzview.setStyle({'stick': {}})
        xyzview.zoomTo()
        showmol(xyzview, height=400, width=500)

        # 提供保存按钮,将 3D 结构写出为 SDF 文件
        if st.button("保存 3D 结构为 SDF"):
            os.makedirs("./result/Prepare_Ligand", exist_ok=True)

            if uploaded_sdf_name:
                # 如果用户上传了 SDF,就使用该 SDF 名（去 .sdf）并加上 _prepared
                base_name = os.path.splitext(uploaded_sdf_name)[0]
                out_filename = base_name + "_prepared.sdf"
            else:
                # 如果没有上传的 SDF,就使用用户输入的文件名（不含 .sdf 后缀）,再加上 .sdf
                if user_defined_filename:
                    out_filename = user_defined_filename.strip() + ".sdf"
                else:
                    # 如果用户也没有输入任何自定义文件名,可给一个默认值
                    out_filename = "ligand_3d.sdf"

            sdf_path = os.path.join("./result/Prepare_Ligand", out_filename)
            writer = Chem.SDWriter(sdf_path)
            writer.write(mol_3d)
            writer.close()
            st.success(f"已将 3D 结构保存到 {sdf_path}")

    # 首先解析用户上传的 SDF
    mol_from_sdf = None
    uploaded_sdf_name = None

    if sdf_file is not None:
        uploaded_sdf_name = sdf_file.name  # 记录用户上传的文件名
        try:
            sdf_supplier = Chem.ForwardSDMolSupplier(sdf_file)
            mols = [m for m in sdf_supplier if m is not None]
            if len(mols) > 0:
                mol_from_sdf = mols[0]
            else:
                st.error("无法从 SDF 文件中解析出分子,请检查文件格式或内容")
        except Exception as e:
            st.error(f"读取 SDF 文件出现错误: {e}")

    if mol_from_sdf:
        # 如果成功解析出上传的 SDF,则展示并保存
        process_and_show_mol(mol_from_sdf, uploaded_sdf_name=uploaded_sdf_name)
    else:
        # 如果用户没有上传 SDF 或上传的 SDF 解析失败,则查看 Ketcher 中有没有输入 SMILES
        if smiles_input:
            mol_from_smiles = Chem.MolFromSmiles(smiles_input)
            if mol_from_smiles:
                user_defined_filename = st.text_input("请输入保存时的 SDF 文件名 (不含`.sdf`): ", value="my_mol")
                process_and_show_mol(
                    mol_from_smiles, 
                    uploaded_sdf_name=None, 
                    user_defined_filename=user_defined_filename
                )
            else:
                st.error("SMILES 无效,请重新输入或确认格式")

# ------------------------------------------------------------------------------
# 口袋预测
# ------------------------------------------------------------------------------
elif page == "口袋预测":
    st.title("口袋预测")

    # 让用户选择如何加载蛋白质
    option = st.radio("Select how to load the protein:", ("上传蛋白质", "加载示例文件"))

    # 用于保存用户上传的蛋白文件名称（用于替换 Pocket Name）
    uploaded_pdb_filename = None

    if option == "上传蛋白质":
        try:
            # 用户上传蛋白质（只出现一次,不会再弹二次上传）
            pdb_file = st.file_uploader("请上传蛋白质文件 (.pdb)", type=["pdb"])
            
            if pdb_file is not None:
                # 记下上传的名称
                uploaded_pdb_filename = pdb_file.name

                # 使用临时文件的方式进行口袋预测
                with tempfile.NamedTemporaryFile(suffix=".pdb", delete=False) as tmp:
                    tmp.write(pdb_file.getvalue())
                    tmp.flush()
                    file_path = tmp.name

                # 调用 p2rank (或其他函数) ,读取该临时文件进行预测
                selected = select_pocket_from_local_protein(
                    file_path,
                    p2rank_home='./utils/p2rank_2.5.1/'
                )
                # 预测完成后删除该临时文件
                os.remove(file_path)

                if selected:
                    pocket = selected
                    st.write("预测到的口袋信息: ", pocket)

                    # 如果 rank=1 的口袋
                    if pocket['rank'] == '1':
                        # 如果上传了文件名,则用之,否则用 pocket['name']
                        final_name = uploaded_pdb_filename if uploaded_pdb_filename else pocket['name']
                        data = {
                            'Pocket Name': [final_name],
                            'Center': [pocket['center']],
                        }
                        df = pd.DataFrame(data)

                        st.write("最优口袋信息预览: ")
                        st.dataframe(df)

                        # 用户点击按钮后,才将CSV保存到指定文件夹
                        if st.button("保存 best_pocket.csv"):
                            os.makedirs("./result/Predict_Pocket", exist_ok=True)
                            csv_path = "./result/Predict_Pocket/best_pocket.csv"
                            df.to_csv(csv_path, index=False)
                            st.success(f"best_pocket.csv 已保存到 {csv_path}")

        except Exception as e:
            st.warning(f"处理上传蛋白时发生错误: {e}")

    elif option == "加载示例文件":
        try:
            # 用示例文件名
            uploaded_pdb_filename = "protein_example.pdb"
            # 调用 p2rank 做预测
            selected = select_pocket_from_local_protein(
                "examples/pocket/protein.pdb", 
                p2rank_home='./utils/p2rank_2.5.1/'
            )
            if selected:
                pocket = selected
                st.write("预测到的口袋信息: ", pocket)

                if pocket['rank'] == '1':
                    data = {
                        'Pocket Name': [uploaded_pdb_filename],
                        'Center': [pocket['center']],
                    }
                    df = pd.DataFrame(data)

                    st.write("最优口袋信息预览: ")
                    st.dataframe(df)

                    if st.button("保存 best_pocket.csv"):
                        os.makedirs("./result/Predict_Pocket", exist_ok=True)
                        csv_path = "./result/Predict_Pocket/best_pocket.csv"
                        df.to_csv(csv_path, index=False)
                        st.success(f"best_pocket.csv 已保存到 {csv_path}")
                        
        except Exception as e:
            st.warning(f"加载示例文件时发生错误: {e}")

# ------------------------------------------------------------------------------
# 分子对接
# ------------------------------------------------------------------------------
elif page == "分子对接":
    import tempfile
    import os
    # import shutil

    st.title("分子对接")
    st.write("请上传蛋白质 (PDB 格式) 和配体 (SDF 格式),并设置对接参数")

    # 让用户上传蛋白质和配体文件
    protein_file = st.file_uploader("上传蛋白质文件 (.pdb)", type=["pdb"])
    ligand_file = st.file_uploader("上传配体文件 (.sdf)", type=["sdf"])

    # 让用户上传口袋预测结果文件（可选）
    st.write("可选：上传口袋预测结果 CSV 文件,将自动填充对接网格参数")
    pocket_csv_file = st.file_uploader("上传口袋预测结果文件 (CSV)", type=["csv"])

    # 初始化 session_state 中的网格参数（如果不存在）
    if 'docking_center_x' not in st.session_state:
        st.session_state.docking_center_x = None
    if 'docking_center_y' not in st.session_state:
        st.session_state.docking_center_y = None
    if 'docking_center_z' not in st.session_state:
        st.session_state.docking_center_z = None
    if 'auto_calculated' not in st.session_state:
        st.session_state.auto_calculated = False
    if 'last_protein_name' not in st.session_state:
        st.session_state.last_protein_name = None

    # 检测蛋白质文件是否改变
    current_protein_name = protein_file.name if protein_file is not None else None
    protein_changed = (current_protein_name != st.session_state.last_protein_name)
    
    if protein_changed and protein_file is not None:
        st.session_state.last_protein_name = current_protein_name
        # 蛋白质文件改变了，重置坐标（如果没有口袋预测）
        if pocket_csv_file is None:
            st.session_state.docking_center_x = None
            st.session_state.docking_center_y = None
            st.session_state.docking_center_z = None

    # 默认网格参数
    center_x = st.session_state.docking_center_x
    center_y = st.session_state.docking_center_y
    center_z = st.session_state.docking_center_z

    # 如果上传了口袋预测结果，优先使用
    if pocket_csv_file is not None:
        try:
            # 读取 CSV 文件并获取中心坐标
            pocket_df = pd.read_csv(pocket_csv_file)
            if "Center" in pocket_df.columns:
                center_coords = pocket_df.loc[0, "Center"]  # 获取第一个口袋的中心坐标
                if isinstance(center_coords, str):
                    coords = [float(c) for c in re.findall(r"[-+]?[0-9]*\.?[0-9]+", center_coords)]
                    if len(coords) == 3:
                        center_x, center_y, center_z = coords
                        st.session_state.docking_center_x = center_x
                        st.session_state.docking_center_y = center_y
                        st.session_state.docking_center_z = center_z
                        st.session_state.auto_calculated = False
                        st.success(f"✅ 已从 CSV 文件读取口袋中心: ({center_x:.2f}, {center_y:.2f}, {center_z:.2f})")
                    else:
                        st.warning("CSV 文件中的 Center 格式不正确,无法自动填充网格参数")
                else:
                    st.warning("CSV 文件中的 Center 格式不正确,无法自动填充网格参数")
            else:
                st.warning("CSV 文件中未找到 Center 列,无法自动填充网格参数")
        except Exception as e:
            st.error(f"读取 CSV 文件时出现错误: {e}")
    
    # 如果用户上传了蛋白质但没有口袋预测结果，自动计算蛋白质质心
    elif protein_file is not None and center_x is None:
        try:
            from biopandas.pdb import PandasPdb
            import io
            
            # 读取蛋白质文件并计算质心
            pdb_content = io.BytesIO(protein_file.getvalue())
            pmol = PandasPdb().read_pdb(pdb_content)
            pdf = pmol.df['ATOM']
            
            if len(pdf) == 0:
                raise ValueError("PDB 文件中没有 ATOM 记录")
            
            center_x = float(pdf['x_coord'].mean())
            center_y = float(pdf['y_coord'].mean())
            center_z = float(pdf['z_coord'].mean())
            
            # 保存到 session_state
            st.session_state.docking_center_x = center_x
            st.session_state.docking_center_y = center_y
            st.session_state.docking_center_z = center_z
            st.session_state.auto_calculated = True
            
            st.info(f"✅ 已自动计算蛋白质质心作为对接网格中心: ({center_x:.2f}, {center_y:.2f}, {center_z:.2f})")
            st.warning("⚠️ 建议先进行'口袋预测'以获得更准确的对接位点！")
        except Exception as e:
            st.error(f"❌ 无法自动计算蛋白质质心: {e}")
            st.warning("请手动输入对接网格参数或先进行口袋预测")
            center_x = 0.0
            center_y = 0.0
            center_z = 0.0
            st.session_state.docking_center_x = center_x
            st.session_state.docking_center_y = center_y
            st.session_state.docking_center_z = center_z
    
    # 设置默认值（如果仍然为 None）
    if center_x is None:
        center_x = 0.0
        center_y = 0.0
        center_z = 0.0
        st.session_state.docking_center_x = center_x
        st.session_state.docking_center_y = center_y
        st.session_state.docking_center_z = center_z

    # 显示网格参数输入框,无论是否上传 CSV 文件
    st.subheader("设置对接口袋参数")
    
    # 添加一个"重置为质心"按钮
    col1, col2 = st.columns([3, 1])
    with col2:
        if protein_file is not None:
            if st.button("🔄 重新计算质心"):
                try:
                    from biopandas.pdb import PandasPdb
                    import io
                    
                    pdb_content = io.BytesIO(protein_file.getvalue())
                    pmol = PandasPdb().read_pdb(pdb_content)
                    pdf = pmol.df['ATOM']
                    
                    center_x = float(pdf['x_coord'].mean())
                    center_y = float(pdf['y_coord'].mean())
                    center_z = float(pdf['z_coord'].mean())
                    
                    st.session_state.docking_center_x = center_x
                    st.session_state.docking_center_y = center_y
                    st.session_state.docking_center_z = center_z
                    st.session_state.auto_calculated = True
                    
                    st.success(f"✅ 已重新计算: ({center_x:.2f}, {center_y:.2f}, {center_z:.2f})")
                except Exception as e:
                    st.error(f"计算失败: {e}")
    
    center_x = st.number_input("Center X", value=float(center_x), format="%.2f", key="input_center_x")
    center_y = st.number_input("Center Y", value=float(center_y), format="%.2f", key="input_center_y")
    center_z = st.number_input("Center Z", value=float(center_z), format="%.2f", key="input_center_z")

    size_x = st.number_input("Size X", value=100.0)
    size_y = st.number_input("Size Y", value=100.0)
    size_z = st.number_input("Size Z", value=100.0)

    # 当用户点击"开始分子对接"时,生成 docking_grid.json 文件并调用对接命令
    if st.button("开始分子对接"):
        # 如果没有上传蛋白质或配体,提示错误
        if not protein_file or not ligand_file:
            st.error("请先上传蛋白质 (pdb) 和配体 (sdf) 文件")
        else:
            # 验证对接网格参数
            if center_x == 0.0 and center_y == 0.0 and center_z == 0.0:
                st.error("⚠️ 警告：对接网格中心为 (0, 0, 0)！")
                st.warning("""
                这可能导致对接失败。请：
                1. 点击上方的 🔄 "重新计算质心" 按钮
                2. 或先进行"口袋预测"以获得准确的对接位点
                3. 或手动输入正确的对接网格中心坐标
                """)
            else:
                st.info(f"📍 对接网格中心: ({center_x:.2f}, {center_y:.2f}, {center_z:.2f})")
            
            try:
                # 创建临时文件夹保存缓存文件
                with tempfile.TemporaryDirectory() as temp_dir:
                    docking_grid = {
                        "center_x": center_x,
                        "center_y": center_y,
                        "center_z": center_z,
                        "size_x": size_x,
                        "size_y": size_y,
                        "size_z": size_z
                    }
                    
                    # 调试信息：打印 docking_grid
                    st.write("调试信息 - 对接网格参数：", docking_grid)
                    docking_grid_path = os.path.join(temp_dir, "docking_grid.json")

                    with open(docking_grid_path, "w") as f:
                        json.dump(docking_grid, f, indent=4)

                    # 保存蛋白质和配体文件到临时目录
                    protein_path = os.path.join(temp_dir, "protein.pdb")
                    ligand_path = os.path.join(temp_dir, "ligand.sdf")

                    with open(protein_path, "wb") as f:
                        f.write(protein_file.getvalue())

                    with open(ligand_path, "wb") as f:
                        f.write(ligand_file.getvalue())

                    # 设置结果保存目录
                    result_dir = "./result/Docking_Result"
                    os.makedirs(result_dir, exist_ok=True)

                    # 构造命令
                    command = (
                        f"conda run -n flashdock python ./models/Uni-Mol/unimol_docking_v2/interface/demo.py "
                        f"--mode single "
                        f"--conf-size 10 "
                        f"--cluster "
                        f"--input-protein {protein_path} "
                        f"--input-ligand {ligand_path} "
                        f"--input-docking-grid {docking_grid_path} "
                        f"--output-ligand-name ligand_predict "
                        f"--output-ligand-dir {result_dir} "
                        f"--steric-clash-fix "
                        f"--model-dir ./models/Uni-Mol/unimol_docking_v2/unimol_docking_v2_240517.pt"
                    )

                    # 执行命令
                    result = subprocess.run(command, shell=True, capture_output=True, text=True)

                    # 根据命令返回值判断是否执行成功
                    if result.returncode == 0:
                        st.success("分子对接完成!")
                        st.text_area("对接输出日志", value=result.stdout, height=150)

                        # 分子对接完成后,处理结果文件
                        try:
                            ligand_output_path = os.path.join(result_dir, "ligand_predict.sdf")

                            # 删除结果目录中除 ligand_predict.sdf 外的所有文件
                            for file_name in os.listdir(result_dir):
                                file_path = os.path.join(result_dir, file_name)
                                if file_name != "ligand_predict.sdf" and os.path.isfile(file_path):
                                    os.remove(file_path)

                            # 重命名 ligand_predict.sdf
                            output_name = f"{os.path.splitext(ligand_file.name)[0]}_{os.path.splitext(protein_file.name)[0]}_docked.sdf"
                            renamed_path = os.path.join(result_dir, output_name)
                            os.rename(ligand_output_path, renamed_path)

                            # 提示用户结果保存位置
                            st.success(f"对接结果保存为 {renamed_path}")

                            # 可视化对接结果
                            st_molstar_docking(
                                protein_path, 
                                renamed_path, 
                                key="5", 
                                height=600
                            )
                        except Exception as e:
                            st.error(f"处理结果文件时出错: {e}")

                    else:
                        st.error("分子对接失败！")
                        
                        # 检查是否是 pocket 为空的错误
                        if "No protein atoms found in the docking grid box" in result.stderr or "Mean of empty slice" in result.stderr:
                            st.error("❌ 在指定的对接网格范围内没有找到蛋白质原子！")
                            st.warning("💡 可能的解决方案：")
                            st.markdown("""
                            1. **建议**：先进行"口袋预测"以获得准确的对接位点
                            2. 检查对接网格的中心坐标是否正确
                            3. 增大对接网格的尺寸（Size X/Y/Z）
                            4. 确认上传的蛋白质文件格式正确
                            """)
                        
                        st.text_area("错误信息", value=result.stderr, height=150)

            except Exception as e:
                st.error(f"对接过程出现错误: {e}")
                st.warning("💡 提示：建议先进行'口袋预测'以获得准确的对接位点！")

# ------------------------------------------------------------------------------
# 批量口袋预测与对接
# ------------------------------------------------------------------------------
elif page == "批量口袋预测与对接":
    import os
    import pandas as pd
    import subprocess
    import tempfile
    import json
    from pathlib import Path
    import streamlit as st
    from streamlit_molstar.pocket import select_pocket_from_local_protein

    st.title("批量口袋预测与分子对接")

    # 定义固定路径
    batch_docking_dir = Path("./batch_docking/input")
    result_dir = Path("./batch_docking/input")
    result_dir.mkdir(parents=True, exist_ok=True)

    # 检查 Batch_Docking 目录是否存在
    if not batch_docking_dir.exists():
        st.error(f"目录 {batch_docking_dir} 不存在请创建该目录并添加 PDB 和 SDF 文件")
    else:
        # 自动生成任务 CSV 文件
        def generate_task_csv():
            pdb_files = list(batch_docking_dir.glob("*.pdb"))
            sdf_files = list(batch_docking_dir.glob("*.sdf"))

            if not pdb_files:
                st.error("在 ./batch_docking/input 文件夹中未找到 PDB 文件请添加至少一个 PDB 文件")
                return None
            if not sdf_files:
                st.error("在 ./batch_docking/input 文件夹中未找到 SDF 文件请添加至少一个 SDF 文件")
                return None

            tasks = []
            for pdb_file in pdb_files:
                for sdf_file in sdf_files:
                    tasks.append({
                        "Protein": pdb_file.name,
                        "Ligand": sdf_file.name,
                        "Run": "Yes"  # 默认所有任务为 "Yes"
                    })

            task_df = pd.DataFrame(tasks)
            return task_df

        # 生成任务 DataFrame
        task_df = generate_task_csv()

        if task_df is not None:
            # 提供下载任务 CSV 的按钮
            csv = task_df.to_csv(index=False).encode('utf-8')
            st.download_button(
                label="下载任务 CSV 文件",
                data=csv,
                file_name="docking_tasks.csv",
                mime="text/csv"
            )

            st.markdown("---")
            st.info("""
                1. 下载上方的任务 CSV 文件
                2. 在本地编辑 CSV 文件,修改 `Run` 列为 `Yes` 的任务将被执行,`No` 列的任务将被跳过
                3. 修改完成后,上传修改后的 CSV 文件并点击“开始批量预测和对接”按钮
            """)

            # 上传修改后的任务 CSV 文件
            uploaded_csv = st.file_uploader("上传修改后的任务 CSV 文件", type=["csv"], key="upload_task_csv")

            if uploaded_csv is not None:
                try:
                    uploaded_tasks_df = pd.read_csv(uploaded_csv)

                    # 检查必要的列是否存在
                    required_columns = {"Protein", "Ligand", "Run"}
                    if not required_columns.issubset(uploaded_tasks_df.columns):
                        st.error(f"上传的任务文件缺少必要的列：{required_columns - set(uploaded_tasks_df.columns)}")
                    else:
                        # 过滤需要运行的任务
                        tasks_to_run = uploaded_tasks_df[uploaded_tasks_df["Run"].str.lower() == "yes"]

                        if tasks_to_run.empty:
                            st.warning("没有任务需要运行,请确保至少有一项任务的 `Run` 列为 `Yes`")
                        else:
                            st.write(f"发现 {len(tasks_to_run)} 个任务需要运行")

                            # 显示需要运行的任务表格
                            st.subheader("待运行的任务列表")
                            st.dataframe(tasks_to_run[['Protein', 'Ligand']].reset_index(drop=True))

                            # 开始批量预测和对接按钮
                            if st.button("开始批量预测和对接", key="start_batch_processing"):
                                log_messages = []
                                progress_bar = st.progress(0)
                                status_text = st.empty()

                                for i, task in tasks_to_run.iterrows():
                                    protein_path = batch_docking_dir / task["Protein"]
                                    ligand_path = batch_docking_dir / task["Ligand"]

                                    # 口袋预测前更新状态
                                    status_text.text(f"任务 {i + 1}/{len(tasks_to_run)}: 正在为 {task['Protein']} 预测口袋...")

                                    # 口袋预测
                                    try:
                                        # 为每个任务传递唯一的 key
                                        pocket_result = select_pocket_from_local_protein(
                                            str(protein_path), 
                                            p2rank_home='./others/p2rank_2.5.1/',
                                            key=f"select_pocket_{i}"  # 添加唯一 key
                                        )
                                    except Exception as e:
                                        log_messages.append(f"任务 {task['Protein']} 和 {task['Ligand']} 口袋预测失败：{e}")
                                        progress_bar.progress((i + 1) / len(tasks_to_run))
                                        continue

                                    if pocket_result:
                                        center_coords = [float(coord) for coord in pocket_result["center"]]
                                        docking_grid = {
                                            "center_x": center_coords[0],
                                            "center_y": center_coords[1],
                                            "center_z": center_coords[2],
                                            "size_x": 100.0,
                                            "size_y": 100.0,
                                            "size_z": 100.0,
                                        }

                                        # 创建临时目录存储对接网格
                                        with tempfile.TemporaryDirectory() as temp_dir:
                                            docking_grid_path = Path(temp_dir) / "docking_grid.json"
                                            with open(docking_grid_path, "w") as f:
                                                json.dump(docking_grid, f, indent=4)

                                            # 更新状态：开始对接
                                            status_text.text(f"任务 {i + 1}/{len(tasks_to_run)}: 正在对接 {task['Protein']} 和 {task['Ligand']}...")

                                            # 构造对接命令
                                            command = (
                                                f"conda run -n flashdock python ./models/Uni-Mol/unimol_docking_v2/interface/demo.py "
                                                f"--mode single "
                                                f"--conf-size 10 "
                                                f"--cluster "
                                                f"--input-protein {protein_path} "
                                                f"--input-ligand {ligand_path} "
                                                f"--input-docking-grid {docking_grid_path} "
                                                f"--output-ligand-name ligand_predict "
                                                f"--output-ligand-dir {result_dir} "
                                                f"--steric-clash-fix "
                                                f"--model-dir ./models/Uni-Mol/unimol_docking_v2/unimol_docking_v2_240517.pt"
                                            )

                                            # 执行对接命令
                                            result = subprocess.run(command, shell=True, capture_output=True, text=True)

                                            if result.returncode == 0:
                                                ligand_output_path = result_dir / "ligand_predict.sdf"
                                                output_name = f"{protein_path.stem}_{ligand_path.stem}_docked.sdf"
                                                renamed_path = result_dir / output_name

                                                try:
                                                    os.rename(ligand_output_path, renamed_path)
                                                    log_messages.append(f"任务 {task['Protein']} 和 {task['Ligand']} 对接完成结果保存为 {renamed_path}")
                                                except Exception as e:
                                                    log_messages.append(f"任务 {task['Protein']} 和 {task['Ligand']} 结果保存失败：{e}")
                                            else:
                                                log_messages.append(f"任务 {task['Protein']} 和 {task['Ligand']} 对接失败错误信息：{result.stderr}")

                                    else:
                                        log_messages.append(f"任务 {task['Protein']} 的口袋信息未找到")

                                    # 更新进度条
                                    progress_bar.progress((i + 1) / len(tasks_to_run))

                                # 所有任务完成后更新状态
                                status_text.text("所有任务已完成")

                                # 显示日志
                                st.success("所有任务已完成")
                                st.text_area("任务日志", value="\n".join(log_messages), height=300)
                except pd.errors.EmptyDataError:
                    st.error("上传的 CSV 文件为空,请检查文件内容")
                except pd.errors.ParserError:
                    st.error("上传的 CSV 文件格式错误,请确保文件为有效的 CSV 格式")
                except Exception as e:
                    st.error(f"读取任务文件时出错：{e}")

# ------------------------------------------------------------------------------
# 预测亲和力
# ------------------------------------------------------------------------------

elif page == "预测亲和力":
    import os
    import tempfile
    import subprocess
    import pandas as pd
    import streamlit as st
    import time
    import matplotlib.pyplot as plt
    import seaborn as sns

    if page == "预测亲和力":
        st.title("预测亲和力")
        st.write("在此页面,你可以进行小分子与蛋白质的结合亲和力预测选择单个预测或批量预测模式")

        # 模式选择
        mode = st.radio("选择模式", ("单个预测", "批量预测"))

        if mode == "单个预测":
            st.subheader("单个蛋白与小分子的亲和力预测")

            # 用户上传蛋白质 PDB 文件
            protein_file = st.file_uploader("上传蛋白质 PDB 文件", type=["pdb"])

            # 用户上传小分子 SDF 文件
            ligand_file = st.file_uploader("上传小分子 SDF 文件", type=["sdf"])

            # 按钮触发预测
            if st.button("开始预测"):
                if protein_file is None:
                    st.error("请上传蛋白质 PDB 文件")
                elif ligand_file is None:
                    st.error("请上传小分子 SDF 文件")
                else:
                    with st.spinner("正在进行亲和力预测,请稍候..."):
                        try:
                            # 创建临时目录
                            with tempfile.TemporaryDirectory() as tmpdir:
                                # 保存上传的蛋白质文件
                                protein_path = os.path.join(tmpdir, protein_file.name)
                                with open(protein_path, "wb") as f:
                                    f.write(protein_file.getbuffer())

                                # 保存上传的小分子文件
                                ligand_path = os.path.join(tmpdir, ligand_file.name)
                                with open(ligand_path, "wb") as f:
                                    f.write(ligand_file.getbuffer())

                                # 输出 CSV 文件路径
                                output_csv_path = os.path.join(tmpdir, "single_prediction.csv")

                                # 调用预测脚本
                                pred_dir = "./models/PLANET"
                                pred_script = "pred.py"
                                pred_script_path = os.path.join(pred_dir, pred_script)

                                cmd = [
                                    "conda run -n flashdock python",
                                    pred_script_path,
                                    "-p", protein_path,
                                    "-l", ligand_path,
                                    "-m", ligand_path,
                                    "-o", output_csv_path
                                ]

                                result = subprocess.run(cmd, capture_output=True, text=True)

                                if result.returncode != 0:
                                    st.error(f"预测过程中发生错误:\n{result.stderr}")
                                else:
                                    if os.path.exists(output_csv_path):
                                        df = pd.read_csv(output_csv_path)
                                        st.success("预测完成! 结果如下: ")
                                        st.dataframe(df)
                                    else:
                                        st.error("预测完成但未找到输出 CSV 文件")
                        except Exception as e:
                            st.error(f"发生异常: {e}")

        elif mode == "批量预测":
            st.subheader("批量蛋白与小分子亲和力预测")

            # 按钮触发预测
            if st.button("开始批量预测"):
                with st.spinner("正在进行批量亲和力预测,请稍候..."):
                    try:
                        batch_dir = "./batch_docking/output"
                        if not os.path.exists(batch_dir):
                            st.error("批量预测目录不存在")
                        else:
                            final_results = []

                            # 扫描文件夹中的 SDF 和 PDB 文件
                            sdf_files = [f for f in os.listdir(batch_dir) if f.endswith(".sdf")]
                            pdb_files = [f for f in os.listdir(batch_dir) if f.endswith(".pdb")]

                            st.write("发现以下蛋白质文件：")
                            st.write(pdb_files)
                            st.write("发现以下配体文件：")
                            st.write(sdf_files)

                            progress_bar = st.progress(0)
                            total_files = len(sdf_files)

                            for i, sdf_file in enumerate(sdf_files):
                                receptor_name = sdf_file.split("_")[0]
                                ligand_name = sdf_file.split("_")[1]
                                pdb_file = os.path.join(batch_dir, receptor_name + ".pdb")
                                sdf_file_path = os.path.join(batch_dir, sdf_file)

                                if os.path.exists(pdb_file):
                                    st.text(f"正在计算第 {i + 1}/{total_files} 对：蛋白 {pdb_file} 和 配体 {sdf_file} 的亲和力...")
                                    with tempfile.TemporaryDirectory() as tmpdir:
                                        output_csv_path_tmp = os.path.join(tmpdir, "temp_result.csv")

                                        cmd = [
                                            "conda run -n flashdock python",
                                            "./models/PLANET/pred.py",
                                            "-p", pdb_file,
                                            "-l", sdf_file_path,
                                            "-m", sdf_file_path,
                                            "-o", output_csv_path_tmp
                                        ]

                                        result = subprocess.run(cmd, capture_output=True, text=True)

                                        if result.returncode == 0 and os.path.exists(output_csv_path_tmp):
                                            temp_df = pd.read_csv(output_csv_path_tmp)
                                            if "Binding_Affinity" in temp_df.columns:
                                                binding_affinity = temp_df["Binding_Affinity"].iloc[0]
                                                final_results.append({
                                                    "Protein_File": receptor_name,
                                                    "Ligand_File": ligand_name,
                                                    "Binding_Affinity": binding_affinity
                                                })
                                        else:
                                            st.error(f"文件 {sdf_file} 处理失败")

                                # 更新进度条
                                progress_bar.progress((i + 1) / total_files)
                                time.sleep(0.1)  # 模拟计算时间

                            if final_results:
                                results_df = pd.DataFrame(final_results)

                                # 保存结果到 Binding_Affinity 文件夹
                                binding_affinity_dir = "./Result/Binding_Affinity"
                                os.makedirs(binding_affinity_dir, exist_ok=True)

                                output_csv_path = os.path.join(binding_affinity_dir, "batch_prediction_results.csv")
                                results_df.to_csv(output_csv_path, index=False)

                                st.success("批量预测完成！结果已保存到以下目录：")
                                st.write(output_csv_path)

                                # 绘制热图
                                st.subheader("亲和力热图")
                                heatmap_data = results_df.pivot(index="Protein_File", columns="Ligand_File", values="Binding_Affinity")
                                plt.figure(figsize=(10, 8), dpi=600)
                                sns.heatmap(heatmap_data, annot=True, cmap="coolwarm", fmt=".2f")
                                plt.xlabel("Ligands")
                                plt.ylabel("Proteins")
                                plt.title("Binding Affinity Heatmap")
                                st.pyplot(plt)

                                heatmap_path = os.path.join(binding_affinity_dir, "binding_affinity_heatmap.png")
                                plt.savefig(heatmap_path, dpi=600)

                                st.write("热图已保存到以下目录：")
                                st.write(heatmap_path)

                            else:
                                st.error("未生成任何预测结果")
                    except Exception as e:
                        st.error(f"发生异常: {e}")




# ------------------------------------------------------------------------------
# 作者信息
# ------------------------------------------------------------------------------
