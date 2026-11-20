"""Streamlit CSAT Calculator V2. Run: streamlit run CSAT_Calculator_V2.py

Property defaults transcribed from the supplied CVOC_properties.png (not live RSL).
Known composition: Calculations 5, 3, 2. Unknown composition: Calculation 4 with PURE solubility.
Assumptions: zero sorption and additive water terms apply to both thresholds.
"""
import io
import re
import math
import numpy as np
import pandas as pd

SOURCE = ('Kueper, B. H., & Davies, K. L. (2009). Assessment and Delineation of '
          'DNAPL Source Zones at Hazardous Waste Sites (EPA/600/R-09/119). U.S. Environmental Protection Agency, '
          'Office of Research and Development, National Risk Management Research Laboratory.')
SOURCE_URL = 'https://nepis.epa.gov/Exe/ZyPURL.cgi?Dockey=P1006Y98.TXT'
# CAS, MW g/mol, density g/mL, solubility mg/L, Koc mL/g, dimensionless H, vapor mm Hg at 25 C.
_RAW = [
 ('Tetrachloroethene (PCE)', '127-18-4',165.85,1.625,210.,265.,.752,18.4),
 ('Trichloroethene (TCE)', '79-01-6',131.4,1.460,1280.,94.,.421,72.),
 ('cis-1,2-Dichloroethene (cis-1,2-DCE)', '156-59-2',96.95,1.284,6400.,36.,.166,175.),
 ('Vinyl Chloride (VC)', '75-01-4',62.5,.912,1100.,19.,1.1,2800.),
 ('1,1,1-Trichloroethane (TCA)', '71-55-6',133.4,1.334,1290.,107.,.711,100.),
 ('1,4-dioxane', '123-91-1',88.11,1.033,None,17.,1.96e-4,29.),
]
PROPERTY_KEYS = ['cas_no','mw','density','solubility','koc','henry','vapor_pressure']
ANALYTE_DB = {r[0]:dict(zip(PROPERTY_KEYS,r[1:])) for r in _RAW}
META = ['location','field_id','lab_id','top_depth','bottom_depth','depth_units','matrix']
REQUIRED = ['location','field_id','top_depth','bottom_depth','depth_units','lab_id','cas_no','analyte','result','units','matrix']
DETAIL = ['field_id','Compound','Observed concentration (mg/kg)',
          'Solubility (mg/L)','Single-component threshold [CTS,i] (mg/kg)','Saturation Index','Threshold met or exceeded']

def norm(x):
    return re.sub(r'[^a-z0-9%]', '',str(x).casefold())

CAS_LOOKUP = {p['cas_no']:n for n,p in ANALYTE_DB.items()}

class ValidationError(ValueError):
    pass

def require(condition, message):
    if not condition:
        raise ValidationError(message)

def normalize_cas(value):
    """Accept ordinary hyphens, Unicode dashes, whitespace and digit-only CAS."""
    text = str(value).strip()
    text = re.sub(r'[\s\u2010-\u2015\u2212-]', '', text)
    if not re.fullmatch(r'\d{5,10}', text):
        return ''
    return f'{text[:-3].lstrip("0") or "0"}-{text[-3:-1]}-{text[-1]}'

def identify(analyte, cas):
    # Chemical identity is CAS-only, regardless of the laboratory's name.
    chemical = CAS_LOOKUP.get(normalize_cas(cas))
    if chemical:
        return chemical
    # An unrecognized chemical CAS must not become a soil property by name.
    if normalize_cas(cas):
        return 'other'
    tokens = re.findall(r'[a-z]+|%', str(analyte).casefold())
    has_moisture = any(t.startswith(('moist', 'moisut')) for t in tokens)
    if has_moisture and not any(t in {'holding', 'capacity'} for t in tokens):
        return 'moisture'
    if ('dry' in tokens and any(t.startswith('den') for t in tokens)
            and ('bulk' in tokens or len(tokens)==2) and 'wet' not in tokens):
        return 'density'
    return 'other'

def validate_input(raw):
    """Fail closed; never interpret a reporting limit or ND as a detected result."""
    df = raw.copy()
    df.columns = [str(c).strip().lower() for c in df.columns]
    require(not df.columns.duplicated().any(), 'Duplicate column names after normalization.')
    require(set(REQUIRED).issubset(df.columns), 'Missing columns: '+', '.join(sorted(set(REQUIRED)-set(df.columns))))
    df = df[REQUIRED].fillna('').copy()
    require(len(df)>0, 'The input file is empty.')
    for c in REQUIRED:
        df[c] = df[c].astype(str).str.strip()
    for c in META+['analyte','result','units']:
        require(df[c].ne('').all(), f'{c} contains blank values. Populate sample metadata on every row.')
    for c in ['top_depth','bottom_depth']:
        df[c] = pd.to_numeric(df[c],errors='coerce')
        require(np.isfinite(df[c]).all(), f'{c} must contain finite numbers.')
    require((df.bottom_depth>=df.top_depth).all(), 'bottom_depth must be >= top_depth.')
    matrix_map={'soil':'Soil','ground water':'Groundwater','groundwater':'Groundwater'}
    df['matrix']=df['matrix'].str.casefold().str.replace(r'\s+', ' ', regex=True).map(matrix_map)
    require(df['matrix'].notna().all(), 'matrix must be Soil, Ground Water, or Groundwater (case-insensitive).')
    df['kind'] = [identify(a,c) for a,c in zip(df.analyte,df.cas_no)]
    warnings = []
    other = df.kind.eq('other')
    groundwater=df.matrix.eq('Groundwater')
    require(df.loc[groundwater,'units'].map(norm).eq('mgl').all(), 'All groundwater concentration records require mg/L units; convert other units before uploading.')
    require(not (groundwater & df.kind.isin(['density','moisture'])).any(), 'Dry bulk density and moisture rows must have matrix Soil.')
    if other.any():
        warnings.append('Unrecognized analytes excluded: '+', '.join(sorted(df.loc[other,'analyte'].unique())))
    relevant = df.loc[~other].copy()
    relevant['value'] = pd.to_numeric(relevant.result,errors='coerce')
    bad = ~np.isfinite(relevant.value) | relevant.value.lt(0)
    require(not bad.any(), 'Results must be finite, nonnegative numeric values (no ND, <, U, or qualifiers). Check rows: '+', '.join(str(i+2) for i in relevant.index[bad]))
    for i,r in relevant.iterrows():
        u = norm(r.units.replace('³','3').replace('^',''))
        valid = u == '%' if r.kind=='moisture' else u in {'gcm3','gcc','gml'} if r.kind=='density' else u==('mgl' if r.matrix=='Groundwater' else 'mgkg')
        require(valid, f'Row {i+2}: {r.analyte} has units {r.units!r}; expected '+('%.' if r.kind=='moisture' else 'g/cm3.' if r.kind=='density' else 'mg/L.' if r.matrix=='Groundwater' else 'mg/kg.'))
    require(not relevant.duplicated(META+['kind']).any(),
            'Duplicate sample/analyte results found. Resolve replicates explicitly; results are not averaged or summed.')
    return df, relevant, warnings

def effective_properties(selected, properties, composition, unknown_mw):
    require(1<=len(selected)<=5 and len(set(selected))==len(selected), 'Select 1–5 unique components.')
    total = sum(composition[n] for n in selected)
    require(all(math.isfinite(composition[n]) and 0<=composition[n]<=100 for n in selected) and total<=100+1e-9,
            'DNAPL composition must be between 0 and 100% in total.')
    require(math.isfinite(unknown_mw) and unknown_mw>0, 'Unknown-component molecular weight must be positive.')
    for n in selected:
        p=properties[n]
        for k in ['mw','solubility','henry']:
            require(p[k] is not None and math.isfinite(p[k]) and p[k]>0, f'{n}: enter a finite positive {k}.')
    unknown = max(0.,100-total)/100
    denom = sum(composition[n]/100/properties[n]['mw'] for n in selected)+unknown/unknown_mw
    require(denom>0,'Mole-fraction denominator is zero.')
    rows=[]
    for n in selected:
        p=properties[n]; mole=composition[n]/100/p['mw']/denom
        rows.append({'Compound':n,'Mole fraction':mole,'Solubility (mg/L)':p['solubility'],
                     'Effective solubility (mg/L)':mole*p['solubility']})
    return pd.DataFrame(rows)

def calculate(raw, selected, properties, effective=None, particle_density=2.65, water_density=1., moisture_basis='Dry mass', mode='unknown'):
    require(all(math.isfinite(x) and x>0 for x in [particle_density,water_density]),'Particle/water densities must be positive.')
    all_rows, data, warnings = validate_input(raw)
    all_rows=all_rows.loc[all_rows.matrix.eq('Soil')]
    data=data.loc[data.matrix.eq('Soil')]
    require(mode in {'known','unknown'}, 'Unknown calculation mode.')
    require(1<=len(selected)<=5 and len(set(selected))==len(selected), 'Select 1–5 unique components.')
    for n in selected:
        for k in ['solubility','henry']:
            value=properties[n][k]
            require(value is not None and math.isfinite(value) and value>0, f'{n}: enter a finite positive {k}.')
    require(mode!='known' or effective is not None, 'Known composition requires effective solubilities.')
    effective_map={} if mode=='unknown' else effective.set_index('Compound')['Effective solubility (mg/L)'].to_dict()
    thresholds=[]; details=[]; summaries=[]; audit=[]
    # Full physical/sample metadata is the key. Never borrow density/moisture across samples.
    for key, group in all_rows.groupby(META,sort=False,dropna=False):
        meta=dict(zip(META,key))
        sub=data.loc[(data[META]==pd.Series(meta)).all(axis=1)]
        label=f"{meta['field_id']} / {meta['lab_id']} / {meta['top_depth']}–{meta['bottom_depth']}"
        vals=sub.set_index('kind')['value'].to_dict()
        missing=[n for n in ['moisture','density']+selected if n not in vals]
        require(not missing, label+': missing '+', '.join(missing)+'. All selected components and soil properties are required.')
        rho=vals['density']; w=vals['moisture']/100
        require(0<rho<=particle_density,label+': dry bulk density must be positive and <= particle density.')
        require(moisture_basis in ['Dry mass','Wet mass'], 'Unknown moisture basis.')
        if moisture_basis=='Wet mass':
            require(w<1,label+': wet-mass moisture must be <100%.')
            w=w/(1-w)
        theta_w=w*rho/water_density
        theta_a=1-rho/particle_density-theta_w
        require(theta_a>1e-12,label+': air-filled porosity is zero or negative; vapor-only thresholds cannot produce a finite CSAT index. Check moisture basis/densities or review the model.')
        ratios=[]
        for n in selected:
            p=properties[n]
            # g/cm3 == kg/L numerically, so (mg/L)/(kg/L) -> mg/kg. No factor 1000.
            vapor_factor=p['henry']*theta_a/rho
            ct=effective_map[n]*vapor_factor if mode=='known' else None
            if mode=='known':
                require(math.isfinite(ct) and ct>=0,label+': invalid effective threshold.')
                comparison = ('Not assessable — zero component fraction' if ct==0 else
                              'Yes' if vals[n]>=ct else 'No')
                threshold_ratio=vals[n]/ct if ct>0 else None
                require(threshold_ratio is None or math.isfinite(threshold_ratio),label+': nonfinite saturation index.')
                thresholds.append({**meta,'Compound':n,'Observed concentration (mg/kg)':vals[n],
                                   'Soil Concentration Threshold [CiT] (mg/kg)':ct,'Saturation Index':threshold_ratio,
                                   'Threshold met or exceeded':comparison})
            if mode=='unknown':
                pure=p['solubility']*vapor_factor
                require(math.isfinite(pure) and pure>0,label+': nonfinite or zero threshold.')
                ratio=vals[n]/pure
                require(math.isfinite(ratio),label+': nonfinite index.')
                ratios.append(ratio)
                details.append({**meta,'Compound':n,'Observed concentration (mg/kg)':vals[n],
                            'Solubility (mg/L)':p['solubility'],
                            'Single-component threshold [CTS,i] (mg/kg)':pure,'Saturation Index':ratio,
                            'Threshold met or exceeded':'Yes' if ratio>=1 else 'No'})
        total=math.fsum(ratios)
        require(math.isfinite(total),label+': CSAT sum is nonfinite.')
        if mode=='unknown':
            summaries.append({**meta,'csat_value':total,'dnapl_present':'Yes' if total>=1 else 'No – attributed to dissolve phase'})
        audit.append({**meta,'dry_bulk_density_g_cm3':rho,'moisture_percent':vals['moisture'],
                      'moisture_basis':moisture_basis,'air_filled_porosity':theta_a})
    if any(n not in selected for n in data.kind if n in ANALYTE_DB):
        warnings.append('Supported analytes outside the selected component list are excluded from results.')
    return (pd.DataFrame(thresholds,columns=THRESHOLD_COLUMNS),pd.DataFrame(details,columns=META+[c for c in DETAIL if c not in META]),
            pd.DataFrame(summaries,columns=META+['csat_value','dnapl_present']),pd.DataFrame(audit),warnings)

GW_RATIO = 'Component ratio [Cobs,i / Si]'
GW_FLAG = 'DNAPL may be present (alpha > 0.01)'
GW_COLUMNS = META+['Compound','Groundwater concentration [Cobs,i] (mg/L)',
                    'Solubility (mg/L)',GW_RATIO,GW_FLAG]
GW_SUMMARY_COLUMNS = META+['alpha',GW_FLAG,'Components included']

def calculate_groundwater(raw, selected, properties):
    """Calculation 6, applied independently to each groundwater sample."""
    all_rows,data,warnings=validate_input(raw)
    require(1<=len(selected)<=5 and len(set(selected))==len(selected), 'Select 1–5 unique components.')
    rows=[]; summaries=[]
    groundwater=all_rows.loc[all_rows.matrix.eq('Groundwater')]
    for key,group in groundwater.groupby(META,sort=False,dropna=False):
        meta=dict(zip(META,key))
        sub=data.loc[(data[META]==pd.Series(meta)).all(axis=1)]
        values=sub.set_index('kind')['value'].to_dict()
        observed=[n for n in selected if n in values]
        label=f"{meta['location']} / {meta['field_id']} / {meta['lab_id']}"
        require(bool(observed),label+': no selected chemical results for Calculation 6.')
        missing=[n for n in selected if n not in values]
        if missing:
            warnings.append(label+': selected components not analyzed: '+', '.join(missing)+'. Alpha uses only measured components and may underestimate DNAPL potential.')
        components=[]
        for name in observed:
            sol=properties[name]['solubility']
            require(sol is not None and math.isfinite(sol) and sol>0, name+': enter a finite positive solubility for Calculation 6.')
            ratio=values[name]/sol
            require(math.isfinite(ratio),label+': nonfinite groundwater component ratio.')
            components.append({'Compound':name,'Groundwater concentration [Cobs,i] (mg/L)':values[name],
                               'Solubility (mg/L)':sol,GW_RATIO:ratio})
        alpha=math.fsum(r[GW_RATIO] for r in components)
        require(math.isfinite(alpha),label+': nonfinite alpha.')
        flag='Yes' if alpha>0.01 else 'No'
        rows.extend({**meta,**r,GW_FLAG:flag} for r in components)
        rows.append({**meta,'Compound':'Sum (alpha)','Groundwater concentration [Cobs,i] (mg/L)':None,
                     'Solubility (mg/L)':None,GW_RATIO:alpha,GW_FLAG:flag})
        summaries.append({**meta,'alpha':alpha,GW_FLAG:flag,'Components included':len(observed)})
    if any(n not in selected for n in data.loc[data.matrix.eq('Groundwater'),'kind'] if n in ANALYTE_DB):
        warnings.append('Supported groundwater analytes outside the selected component list are excluded from alpha.')
    return pd.DataFrame(rows,columns=GW_COLUMNS),pd.DataFrame(summaries,columns=GW_SUMMARY_COLUMNS),warnings

def excel_bytes(sheets):
    out=io.BytesIO()
    with pd.ExcelWriter(out,engine='xlsxwriter',engine_kwargs={'options':{'strings_to_formulas':False,'strings_to_urls':False}}) as writer:
        for name,frame in sheets.items():
            frame.to_excel(writer,sheet_name=name,index=False)
            ws=writer.sheets[name]; ws.freeze_panes(1,0)
            if len(frame.columns):
                ws.autofilter(0,0,len(frame),len(frame.columns)-1)
            for j,c in enumerate(frame.columns):
                ws.set_column(j,j,min(60,max(14,len(c)+2)))
    return out.getvalue()

def read_input(uploaded, sheet=0):
    if uploaded.name.lower().endswith('.csv'):
        return pd.read_csv(uploaded,dtype=str,keep_default_na=False)
    return pd.read_excel(uploaded,sheet_name=sheet,dtype=str,keep_default_na=False)

THRESHOLD_COLUMNS = META+['Compound','Observed concentration (mg/kg)','Soil Concentration Threshold [CiT] (mg/kg)','Saturation Index','Threshold met or exceeded']

def filtered_table(st, frame, columns, prefix, disabled=False):
    a,b=st.columns(2)
    fields=a.multiselect('Filter field sample ID',list(frame.field_id.unique()),key=prefix+'_fields',disabled=disabled)
    names=b.multiselect('Filter analyte',list(frame.Compound.unique()),key=prefix+'_names',disabled=disabled)
    view=frame
    if fields: view=view.loc[view.field_id.isin(fields)]
    if names: view=view.loc[view.Compound.isin(names)]
    st.dataframe(view[columns],hide_index=True,width='stretch')
    st.caption('Empty filters show all rows. Downloads include all samples and analytes, independent of display filters.')

def main():
    import streamlit as st
    st.set_page_config(page_title='DNAPL CSAT Calculator V2',layout='wide')
    st.title('DNAPL CSAT Calculator · Version 2')
    st.markdown('''<style>
    .st-key-threshold_inactive, .st-key-csat_inactive, .st-key-groundwater_inactive {
        opacity:0.48; pointer-events:none;
    }
    </style>''',unsafe_allow_html=True)
    uploaded=st.file_uploader('Upload input file (CSV or Excel)',type=['csv','xlsx','xls'])
    raw=None
    if uploaded is not None:
        try:
            sheet=0
            if not uploaded.name.lower().endswith('.csv'):
                book=pd.ExcelFile(uploaded)
                sheet=st.selectbox('Input worksheet',book.sheet_names)
            uploaded.seek(0)
            raw=read_input(uploaded,sheet)
        except Exception as exc:
            st.error(f'Could not read input file: {exc}')
    st.caption('Required columns: '+', '.join(REQUIRED)+'. Chemical records match by cas_no only. matrix accepts Soil, Ground Water, or Groundwater. Groundwater concentrations must be mg/L; soil chemicals must be mg/kg. Moisture and dry bulk density are separate Soil rows.')
    st.header('Parameters')
    known=st.radio('Is the DNAPL composition by mass known or estimated?', ['Yes','No'],horizontal=True)=='Yes'
    st.caption('Yes: Calculations 5 → 3 → 2, then compare measured concentrations with component thresholds. No: Calculation 4 for soil and Calculation 6 for groundwater, using solubility values. Both workflows remain visible.')
    st.subheader('Analyte Properties')
    count=st.number_input('Number of DNAPL components',min_value=1,max_value=5,value=1,step=1)
    names=list(ANALYTE_DB); selected=[]; properties={}
    for i in range(count):
        n=st.selectbox(f'Component {i+1}',names,index=i,key=f'component_{i}')
        selected.append(n); defaults=ANALYTE_DB[n]; p=defaults.copy()
        with st.container(border=True):
            st.markdown(f'**{n}** · CAS {defaults["cas_no"]}')
            cols=st.columns(3)
            p['mw']=cols[0].number_input('Molecular weight (g/mol)',min_value=.01,value=defaults['mw'],format='%.2f',disabled=not known,key=n+str(i)+'mw')
            p['solubility']=cols[1].number_input('Solubility (mg/L)',min_value=.01,value=defaults['solubility'],format='%.2f',key=n+str(i)+'sol')
            p['henry']=cols[2].number_input('Henry’s constant (dimensionless)',min_value=.000000001,value=defaults['henry'],format='%.9f',key=n+str(i)+'henry')
            if defaults['solubility'] is None:
                st.info('The supplied table lists 1,4-dioxane as miscible. Enter a technically justified numeric solubility; no numeric default is assumed.')
        properties[n]=p
    st.caption('Property defaults: supplied CVOC_properties screenshot. Molecular weight and solubility display two decimal places; Henry’s constant retains its precision. Molecular weight is inactive when composition is unknown.')
    st.subheader('Site Properties')
    signature=tuple(selected)
    if st.session_state.get('_composition_components') != signature:
        for n in dict.fromkeys(selected):
            st.session_state['mass_'+n]=100./count
        st.session_state['_composition_components']=signature
    composition={n:st.number_input(f'DNAPL composition by mass (%) — {n}',min_value=0.,max_value=100.,step=1.,format='%.2f',disabled=not known,key='mass_'+n) for n in dict.fromkeys(selected)}
    total=sum(composition.values())
    st.caption(f'Entered composition: {total:.2f}%; unknown balance: {max(0.,100-total):.2f}%. Composition is ignored in the No workflow. Changing components resets percentages to equal shares.')
    mw_estimate=float(np.mean([properties[n]['mw'] for n in selected]))
    unknown_mw=st.number_input('Estimated molecular weight of unknown balance (g/mol)',min_value=.01,value=round(mw_estimate,2),format='%.2f',disabled=not known or total>=100-1e-9)
    st.caption('For a composition below 100%, the mole-fraction denominator includes the unknown balance using this estimated molecular weight. Known components are not renormalized.')
    a,b,c=st.columns(3)
    particle=a.number_input('Soil particle density (g/cm³)',min_value=.000001,value=2.65,format='%.6f')
    water=b.number_input('Water density (g/cm³)',min_value=.000001,value=1.,format='%.6f')
    basis=c.selectbox('Laboratory moisture reporting basis',['Dry mass','Wet mass'])
    st.caption('Dry mass: water / dry soil. Wet mass: water / wet soil, converted to dry mass for porosity. Dry bulk density is used without a second moisture correction.')
    eff=pd.DataFrame(columns=['Compound','Mole fraction','Solubility (mg/L)','Effective solubility (mg/L)'])
    thresholds=pd.DataFrame(columns=THRESHOLD_COLUMNS)
    detail=pd.DataFrame(columns=META+[c for c in DETAIL if c not in META])
    summary=pd.DataFrame(columns=META+['csat_value','dnapl_present'])
    audit=pd.DataFrame(); ready=False; gw_ready=False; has_groundwater=False
    gw_table=pd.DataFrame(columns=GW_COLUMNS)
    gw_summary=pd.DataFrame(columns=GW_SUMMARY_COLUMNS)
    if raw is not None:
        normalized=raw.rename(columns=lambda c:str(c).strip().lower())
        if 'matrix' in normalized and not normalized.columns.duplicated().any():
            has_groundwater=normalized['matrix'].astype(str).str.strip().str.casefold().str.replace(r'\s+', '', regex=True).eq('groundwater').any()
    valid_input=False
    try:
        require(len(set(selected))==len(selected),'Select unique components.')
        if known:
            eff=effective_properties(selected,properties,composition,unknown_mw)
        if raw is not None:
            validate_input(raw)
            valid_input=True
    except ValidationError as exc:
        st.error(str(exc))
    if valid_input:
        try:
            thresholds,detail,summary,audit,warnings=calculate(raw,selected,properties,eff if known else None,particle,water,basis,mode='known' if known else 'unknown')
            for message in warnings: st.warning(message)
            ready=not (thresholds.empty if known else detail.empty)
        except ValidationError as exc:
            st.error('Soil calculation: '+str(exc))
        if not known and has_groundwater:
            try:
                gw_table,gw_summary,warnings=calculate_groundwater(raw,selected,properties)
                for message in warnings: st.warning(message)
                gw_ready=not gw_table.empty
            except ValidationError as exc:
                st.error('Groundwater calculation: '+str(exc))
    if raw is None:
        st.info('Upload an input file to calculate sample results. Both workflows are shown below.')
    settings=pd.DataFrame([
        {'setting':'composition_mode','value':'Known/estimated' if known else 'Unknown'},
        {'setting':'particle_density_g_cm3','value':particle},
        {'setting':'water_density_g_cm3','value':water},
        {'setting':'moisture_basis','value':basis},
        {'setting':'Kd_rho_b','value':0}, {'setting':'additive_theta_w','value':0},
        {'setting':'unknown_balance_percent','value':max(0.,100-total) if known else 'Not used'},
        {'setting':'unknown_mw_g_mol','value':unknown_mw if known else 'Not used'},
        {'setting':'source','value':SOURCE}, {'setting':'source_url','value':SOURCE_URL}])
    props=pd.DataFrame([{'Compound':n,'cas_no':properties[n]['cas_no'],
        'Molecular weight (g/mol)':properties[n]['mw'],
        'Solubility (mg/L)':properties[n]['solubility'],
        'Henry constant':properties[n]['henry'],
        **({'composition_percent':composition[n]} if known else {})} for n in selected])
    common={'Sample Audit':audit,'Analyte Properties':props,'Settings':settings}
    if valid_input and not ready:
        st.info('No soil results are available for the active soil workflow. Soil calculations do not use groundwater concentrations.')
    with st.container(border=True,key='threshold_active' if known else 'threshold_inactive'):
        st.header('Threshold Chemical Concentration in Soil Based on Partitioning Relationships Where Composition of DNAPL is Known')
        st.subheader('Effective Solubility')
        st.caption('Active — known/estimated composition' if known else 'Inactive — select Yes to use Calculations 5 and 3.')
        st.latex(r'm_i=\frac{ms_i/mw_i}{\sum_j ms_j/mw_j},\qquad C_i=m_iS_i')
        st.dataframe(eff,hide_index=True,width='stretch',column_config={
            'Solubility (mg/L)':st.column_config.NumberColumn(format='%.2f'),
            'Effective solubility (mg/L)':st.column_config.NumberColumn(format='%.2f')})
        st.caption('mi: mole fraction; msi: mass fraction (composition % / 100); mwi: molecular weight (g/mol); j: all known components plus estimated unknown balance; Ci: effective solubility (mg/L); Si: solubility (mg/L).')
        st.caption('Calculations 3 and 5. '+SOURCE)
        st.subheader('Component Soil Thresholds')
        st.caption('Active — compare each component separately; ratios are not summed.' if known else 'Inactive — select Yes for composition-adjusted thresholds.')
        filtered_table(st,thresholds,['field_id','Compound','Observed concentration (mg/kg)','Soil Concentration Threshold [CiT] (mg/kg)','Saturation Index','Threshold met or exceeded'],'threshold',disabled=not known or not ready)
        st.latex(r'\theta_a=1-\rho_b/\rho_s-w\rho_b/\rho_w,\qquad C_i^T=C_iH_i^\prime\theta_a/\rho_b')
        st.caption('Assumptions: Kd × ρb = 0 and additive θw = 0. Measured moisture still reduces air-filled porosity. θa: air-filled porosity; ρb: dry bulk density (g/cm³ = kg/L); ρs: particle density; ρw: water density; w: dry-mass moisture fraction; H′: dimensionless Henry’s constant; Ci: effective solubility (mg/L); CiT: soil threshold (mg/kg).')
        st.caption('Calculation 2, using the stated assumptions. '+SOURCE)
        st.latex(r'\mathrm{Saturation\ Index}_i=C^T_{obs,i}/C_i^T')
        st.caption('Saturation Index = measured soil concentration / composition-adjusted soil threshold. Below 1: below threshold; 1: at threshold; above 1: above threshold. Compare each component separately; do not sum these ratios or interpret their sum as Calculation 4 CSAT.')
        st.caption('Yes indicates a screening exceedance under the stated assumptions, not confirmed DNAPL presence. No does not establish absence. Excluding sorbed and dissolved contributions lowers the threshold and can produce a screening exceedance without DNAPL. A zero component fraction gives a zero threshold: the ratio is left blank and the comparison is marked not assessable; review the assumed composition if that component is measured.')
        st.download_button('Export composition-adjusted thresholds (Excel)',excel_bytes({'Component Thresholds':thresholds,'Effective Solubility':eff,**common}) if known and ready else b'',file_name='CSAT_V2_Known_Composition.xlsx',disabled=not known or not ready)
    with st.container(border=True,key='csat_inactive' if known else 'csat_active'):
        st.header('Threshold Chemical Concentration in Soil Based on Partitioning Relationships Where Composition of DNAPL is Not Known')
        st.caption('Inactive — select No to use Calculation 4.' if known else 'Active — unknown composition')
        st.latex(r'\sum_{i=1}^{n}\frac{C^T_{obs,i}}{C^T_{S,i}}\geq 1')
        st.info('Calculation 4 uses solubility and does not require composition percentages. Assumptions: sorption and additive water terms are zero. Measured moisture still determines air-filled porosity using θa = 1 − ρb/ρs − wρb/ρw.')
        filtered_table(st,detail,DETAIL,'csat',disabled=known or not ready)
        st.caption('CTobs,i: reported soil concentration of component i (mg/kg); CTS,i: single-component soil partitioning threshold of component i (mg/kg), calculated using Calculation 2 with single-component solubility as input; n: number of selected measured components included in the soil sample. The table Saturation Index is the component concentration divided by its single-component threshold. The sample total is the sum of those indices.')
        st.caption('Threshold met or exceeded compares each measured component concentration with its single-component threshold [CTS,i] (component index ≥ 1). The sample CSAT sum can reach 1 even if every individual component is below its threshold.')
        st.caption('Calculation 4. '+SOURCE)
        st.caption('CSAT ≥ 1 gives Yes; CSAT < 1 gives the requested “No – attributed to dissolve phase” export label. These are screening classifications under the modified model; a value below 1 does not establish absence. Numeric chemical results must be detected concentrations on a dry-soil basis.')
        a,b=st.columns(2)
        a.download_button('Export Non-Composition-Adjusted Thresholds (Excel)',excel_bytes({'Non-Adjusted Thresholds':detail[DETAIL],'Sample Metadata':detail,**common}) if not known and ready else b'',file_name='CSAT_V2_Non_Composition_Adjusted_Thresholds.xlsx',disabled=known or not ready)
        b.download_button('Export CSAT sum by sample (Excel)',excel_bytes({'CSAT by Sample':summary,**common}) if not known and ready else b'',file_name='CSAT_V2_By_Sample.xlsx',disabled=known or not ready)

    gw_active=not known and has_groundwater and gw_ready
    gw_settings=pd.DataFrame([
        {'setting':'calculation','value':'6'},
        {'setting':'composition_mode','value':'Unknown'},
        {'setting':'groundwater_concentration_units','value':'mg/L'},
        {'setting':'alpha_rule','value':'> 0.01 (strict)'},
        {'setting':'assumptions','value':'Equal borehole dilution, dispersion and degradation effects for each component'},
        {'setting':'source','value':SOURCE}, {'setting':'source_url','value':SOURCE_URL}])
    with st.container(border=True,key='groundwater_active' if gw_active else 'groundwater_inactive'):
        st.header('1% Effective Solubility Threshold Not Knowing DNAPL Composition')
        if gw_active:
            st.caption('Active — unknown composition and valid groundwater results available.')
        elif known:
            st.caption('Inactive — select No for DNAPL composition and upload groundwater records to activate this section.')
        elif not has_groundwater:
            st.caption('Inactive — upload records with matrix Ground Water or Groundwater to activate this section.')
        else:
            st.caption('Inactive — resolve the input or calculation errors above to activate this section.')
        st.latex(r'\alpha=\sum_{i=1}^{n}\frac{C_i^{obs}}{S_i},\qquad \mathrm{DNAPL\ may\ be\ present\ if}\ \alpha>0.01')
        st.caption('Groundwater only. The Yes/No flag uses the sample alpha (sum), not any single component ratio. Sum rows are included for each sample; No means the 1% rule was not exceeded and does not establish DNAPL absence.')
        options=list(gw_table.field_id.unique())
        fields=st.multiselect('Filter groundwater field sample ID',options,key='gw_fields',disabled=not gw_active)
        view=gw_table.loc[gw_table.field_id.isin(fields)] if fields else gw_table
        st.dataframe(view[['field_id','Compound','Groundwater concentration [Cobs,i] (mg/L)','Solubility (mg/L)',GW_RATIO,GW_FLAG]],hide_index=True,width='stretch')
        st.caption('Empty filter shows all samples. Downloads include all groundwater samples, independent of the filter. Full sample identifiers are preserved in the export; samples are not pooled across field IDs or depths.')
        st.download_button('Export Groundwater 1% Threshold (Excel)',excel_bytes({'Groundwater Components':gw_table,'Groundwater by Sample':gw_summary,'Analyte Properties':props,'Settings':gw_settings}) if gw_active else b'',file_name='CSAT_V2_Groundwater_1_Percent_Threshold.xlsx',disabled=not gw_active)
        st.caption('Calculation 6 evaluates possible DNAPL presence near the sampled monitoring location using the cumulative solubility-normalized groundwater concentrations. Assumptions: borehole dilution, dispersion, and degradation affect each component equally. Each sample is assessed separately; spatially uniform DNAPL composition is not required.')
        st.caption('Cobs,i: sampled groundwater concentration of component i (mg/L); Si: single-component solubility of component i (mg/L); alpha (α): dimensionless cumulative mole-fraction estimate for the sample set, calculated as the sum of Cobs,i / Si; n: number of selected, measured components included in that groundwater sample. The rule uses a strict > 0.01 comparison on unrounded alpha.')
        st.caption('Missing or unmeasured DNAPL components and aqueous-phase degradation can produce false negatives or underestimate DNAPL potential. Use exceedances as screening evidence of possible nearby DNAPL; do not use a non-exceedance to demonstrate absence.')
        st.caption('Calculation #6 in '+SOURCE)

if __name__=='__main__':
    main()
