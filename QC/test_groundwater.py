import copy
import io
import unittest
from pathlib import Path
import numpy as np
import pandas as pd
import CSAT_Calculator_V2 as app

class GroundwaterTests(unittest.TestCase):
    def setUp(self):
        self.props=copy.deepcopy(app.ANALYTE_DB)
        self.names=list(self.props)[:2]
        self.soil=pd.read_csv(Path(__file__).with_name('example_input.csv'),dtype=str,keep_default_na=False)
        self.gw=self.soil.loc[self.soil.cas_no.ne('')].copy()
        self.gw['matrix']='Ground Water'; self.gw['units']='mg/L'
        self.gw.loc[self.gw.cas_no.eq('127-18-4'),'result']='1.89' # .009 of solubility
        self.gw.loc[self.gw.cas_no.eq('79-01-6'),'result']='5.12' # .004 of solubility
    def calc(self,raw=None):
        return app.calculate_groundwater(self.gw if raw is None else raw,self.names,self.props)
    def test_sum_rule_and_layout(self):
        detail,summary,_=self.calc()
        np.testing.assert_allclose(summary.alpha,.013)
        self.assertTrue(summary[app.GW_FLAG].eq('Yes').all())
        self.assertEqual(len(detail),6)
        self.assertEqual(detail.Compound.eq('Sum (alpha)').sum(),2)
        self.assertTrue((detail.loc[detail.Compound.ne('Sum (alpha)'),app.GW_RATIO]<.01).all())
        self.assertEqual(list(detail),app.GW_COLUMNS)
    def test_screenshot_arithmetic(self):
        names=list(self.props)[:5]
        concentrations=[4.4,1.8,3.5,4.,48.]
        solubilities=[1100.,200.,500.,500.,8000.]
        rows=[]
        base=self.gw.iloc[0].to_dict()
        for name,conc,sol in zip(names,concentrations,solubilities):
            self.props[name]['solubility']=sol
            rows.append({**base,'cas_no':self.props[name]['cas_no'],'analyte':'lab name','result':conc})
        detail,summary,_=app.calculate_groundwater(pd.DataFrame(rows),names,self.props)
        self.assertAlmostEqual(summary.iloc[0].alpha,.034)
        np.testing.assert_allclose(detail.iloc[:-1][app.GW_RATIO],[.004,.009,.007,.008,.006])
        self.assertEqual(detail.iloc[-1][app.GW_FLAG],'Yes')
    def test_strict_one_percent_boundary(self):
        self.props[self.names[0]]['solubility']=100.
        self.names=self.names[:1]
        raw=self.gw.loc[self.gw.cas_no.eq('127-18-4')].copy()
        for conc,flag in [('0.999999','No'),('1','No'),('1.000001','Yes')]:
            raw['result']=conc
            self.assertTrue(self.calc(raw)[1][app.GW_FLAG].eq(flag).all())
    def test_matrix_required_and_allowed(self):
        with self.assertRaises(app.ValidationError):self.calc(self.gw.drop(columns='matrix'))
        for spelling in ['Groundwater','ground water','GROUNDWATER','  Ground Water  ']:
            raw=self.gw.copy(); raw['matrix']=spelling
            self.assertTrue(self.calc(raw)[1].matrix.eq('Groundwater').all())
        for spelling in ['water','sediment','','soil water']:
            raw=self.gw.copy();raw['matrix']=spelling
            with self.assertRaises(app.ValidationError):self.calc(raw)
        raw=self.soil.copy();raw['matrix']='sOiL'
        self.assertTrue(app.validate_input(raw)[0].matrix.eq('Soil').all())
    def test_groundwater_units_even_unrecognized(self):
        for units in ['ug/L','mg/kg','%']:
            raw=self.gw.copy();raw.loc[raw.index[0],'units']=units
            with self.assertRaises(app.ValidationError):self.calc(raw)
        raw=self.gw.copy();raw.loc[raw.index[0],'cas_no']='999-99-9';raw.loc[raw.index[0],'units']='ug/L'
        with self.assertRaises(app.ValidationError):self.calc(raw)
        raw=self.gw.copy();raw['units']='MG / L';self.calc(raw)
    def test_mixed_matrices_keep_samples_separate(self):
        raw=pd.concat([self.soil,self.gw],ignore_index=True)
        unknown=app.calculate(raw,self.names,self.props)[2]
        np.testing.assert_allclose(unknown.csat_value,app.calculate(self.soil,self.names,self.props)[2].csat_value)
        self.assertTrue(unknown.matrix.eq('Soil').all())
        np.testing.assert_allclose(self.calc(raw)[1].alpha,self.calc()[1].alpha)
    def test_groundwater_needs_no_soil_properties_or_henry(self):
        for name in self.names:self.props[name]['henry']=None
        self.assertEqual(len(self.calc()[1]),2)
    def test_missing_panel_warns_without_imputed_zero(self):
        raw=self.gw.loc[self.gw.cas_no.eq('127-18-4')]
        detail,summary,warnings=self.calc(raw)
        np.testing.assert_allclose(summary.alpha,.009)
        self.assertTrue(summary['Components included'].eq(1).all())
        self.assertTrue(any('not analyzed' in w for w in warnings))
        self.assertEqual(len(detail),4)
    def test_duplicates_after_matrix_normalization(self):
        duplicate=self.gw.iloc[[0]].copy();duplicate['matrix']='groundwater'
        with self.assertRaises(app.ValidationError):self.calc(pd.concat([self.gw,duplicate],ignore_index=True))
    def test_groundwater_invalid_solubility_and_qualified_results(self):
        for value in [None,0.,float('nan')]:
            self.props[self.names[0]]['solubility']=value
            with self.assertRaises(app.ValidationError):self.calc()
        self.props=copy.deepcopy(app.ANALYTE_DB)
        for result in ['ND','<1','nan','-1']:
            raw=self.gw.copy();raw.loc[raw.index[0],'result']=result
            with self.assertRaises(app.ValidationError):self.calc(raw)
    def test_unknown_soil_flag(self):
        detail=app.calculate(self.soil,self.names,self.props)[1]
        expected=detail['Observed concentration (mg/kg)']>=detail['Single-component threshold [CTS,i] (mg/kg)']
        self.assertEqual(list(detail['Threshold met or exceeded']),['Yes' if x else 'No' for x in expected])
    def test_groundwater_workbook(self):
        detail,summary,_=self.calc()
        workbook=app.excel_bytes({'Groundwater Components':detail,'Groundwater by Sample':summary})
        actual=pd.read_excel(io.BytesIO(workbook),sheet_name='Groundwater by Sample')
        np.testing.assert_allclose(actual.alpha,summary.alpha)
        self.assertEqual(list(actual),app.GW_SUMMARY_COLUMNS)
    def test_soil_only_has_no_groundwater_results(self):
        detail,summary,_=self.calc(self.soil)
        self.assertTrue(detail.empty and summary.empty)

if __name__=='__main__':unittest.main()
