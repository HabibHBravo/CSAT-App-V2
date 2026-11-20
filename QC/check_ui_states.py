"""Streamlit UI integration checks. Run: python check_ui_states.py"""
from pathlib import Path
from streamlit.testing.v1 import AppTest

ROOT=Path(__file__).resolve().parent
GW_TITLE='1% Effective Solubility Threshold Not Knowing DNAPL Composition'

def loaded(name):
    code=f'''
import sys,io
sys.path.insert(0,{str(ROOT)!r})
import streamlit as st
from unittest.mock import patch
import CSAT_Calculator_V2 as app
class Upload(io.BytesIO):
    name='input.csv'
raw=Upload(open({str(ROOT/name)!r},'rb').read())
with patch.object(st,'file_uploader',return_value=raw):
    app.main()
'''
    at=AppTest.from_string(code).run()
    at.number_input[0].set_value(2).run()
    assert not at.exception,at.exception
    return at

def has_groundwater_card(at):
    return any(h.value==GW_TITLE for h in at.header)

# Soil-only: groundwater card stays visible and inactive in both workflows.
soil=loaded('example_input.csv')
assert has_groundwater_card(soil) and soil.multiselect[-1].disabled
assert 'Saturation Index' in soil.dataframe[1].value
assert 'Soil Concentration Threshold [CiT] (mg/kg)' in soil.dataframe[1].value
soil.radio[0].set_value('No').run()
assert not soil.exception and not soil.error and has_groundwater_card(soil) and soil.multiselect[-1].disabled
assert 'Single-component threshold [CTS,i] (mg/kg)' in soil.dataframe[2].value
assert 'Threshold met or exceeded' in soil.dataframe[2].value
assert [b.proto.disabled for b in soil.get('download_button')]==[True,False,False,True]
# Groundwater-only: no soil results or missing-density errors.
gw=loaded('example_groundwater_input.csv')
assert has_groundwater_card(gw) and gw.multiselect[-1].disabled
gw.radio[0].set_value('No').run()
assert not gw.exception and not gw.error
assert has_groundwater_card(gw)
assert len(gw.header)==4 and gw.header[-1].value==GW_TITLE
assert [b.proto.disabled for b in gw.get('download_button')]==[True,True,True,False]
assert len(gw.dataframe[-1].value)==6
assert len(gw.dataframe[2].value)==0
assert 'Sum (alpha)' in list(gw.dataframe[-1].value.Compound)
gw.multiselect[-1].set_value(['MW-01-GW']).run()
assert len(gw.dataframe[-1].value)==3
assert set(gw.dataframe[-1].value['DNAPL may be present (alpha > 0.01)'])=={'Yes'}
gw.radio[0].set_value('Yes').run()
assert not gw.exception and has_groundwater_card(gw) and gw.multiselect[-1].disabled
# Mixed input: soil and groundwater exports can be enabled together.
mixed=loaded('example_mixed_input.csv')
mixed.radio[0].set_value('No').run()
assert not mixed.exception and not mixed.error and has_groundwater_card(mixed)
assert [b.proto.disabled for b in mixed.get('download_button')]==[True,False,False,False]
assert len(mixed.dataframe[2].value)==4 and len(mixed.dataframe[3].value)==6
# Bad units display an error and block exports; the eligible card remains visible.
code=f'''
import sys,io,pandas as pd
sys.path.insert(0,{str(ROOT)!r})
import streamlit as st
from unittest.mock import patch
import CSAT_Calculator_V2 as app
raw=pd.read_csv({str(ROOT/'example_groundwater_input.csv')!r},dtype=str)
raw['units']='ug/L'
class Upload(io.BytesIO):name='bad.csv'
upload=Upload(raw.to_csv(index=False).encode())
with patch.object(st,'file_uploader',return_value=upload):app.main()
'''
bad=AppTest.from_string(code).run()
bad.radio[0].set_value('No').run()
assert not bad.exception and bad.error and has_groundwater_card(bad)
assert all(b.proto.disabled for b in bad.get('download_button'))
# Missing uploads retain all three calculation boxes, with exports disabled.
empty=AppTest.from_file(str(ROOT/'CSAT_Calculator_V2.py')).run()
assert not empty.exception and len(empty.header)==4 and has_groundwater_card(empty) and empty.multiselect[-1].disabled
assert all(b.proto.disabled for b in empty.get('download_button'))
# Check visible UI text does not contain retired terminology.
for at in [soil,gw,mixed,empty]:
    labels=' '.join(str(getattr(w,'label','')) for w in at.number_input)
    texts=' '.join(str(w.value) for w in list(at.caption)+list(at.info))
    assert 'handbook solubility' not in (labels+texts).casefold()
print('UI checks passed: soil, groundwater and mixed inputs, always-visible card activation, filters, exports, labels, invalid units and no upload.')
