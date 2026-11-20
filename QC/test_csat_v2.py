import copy
import io
import unittest
from pathlib import Path
import pandas as pd
import numpy as np
import CSAT_Calculator_V2 as app

class CalculatorTests(unittest.TestCase):
    def setUp(self):
        self.raw=pd.read_csv(Path(__file__).with_name('example_input.csv'),dtype=str,keep_default_na=False)
        self.names=list(app.ANALYTE_DB)[:2]
        self.props=copy.deepcopy(app.ANALYTE_DB)
        self.mass={n:50. for n in self.names}
        self.eff=app.effective_properties(self.names,self.props,self.mass,150.)
    def calc(self,raw=None,eff=None,**kw):
        return app.calculate(self.raw if raw is None else raw,self.names,self.props,self.eff if eff is None else eff,**kw)
    def test_hand_computed_values(self):
        th,dt,su,audit,_=self.calc()
        air=1-1.4/2.65-.2*1.4
        self.assertAlmostEqual(air,.1916981132075472)
        self.assertAlmostEqual(audit.iloc[0].air_filled_porosity,air)
        pce=self.names[0]; tce=self.names[1]
        expected=(155/(210*.752*air/1.4)+145/(1280*.421*air/1.4))
        self.assertAlmostEqual(su.iloc[0].csat_value,expected)
        self.assertEqual(su.iloc[0].dnapl_present,'Yes')
        self.assertEqual(su.iloc[1].dnapl_present,'No – attributed to dissolve phase')
        mole=(.5/165.85)/(.5/165.85+.5/131.4)
        self.assertAlmostEqual(self.eff.iloc[0]['Mole fraction'],mole)
        known=self.calc(mode='known')
        self.assertAlmostEqual(known[0].iloc[0]['Soil Concentration Threshold [CiT] (mg/kg)'],mole*210*.752*air/1.4)
        self.assertTrue(th.empty)
        self.assertTrue(known[1].empty)
        self.assertTrue(known[2].empty)
        self.assertEqual(list(su),app.META+['csat_value','dnapl_present'])
    def test_composition_does_not_change_csat(self):
        eff=app.effective_properties(self.names,self.props,dict(zip(self.names,[10.,90.])),150.)
        a=self.calc(); b=self.calc(eff=eff)
        np.testing.assert_allclose(a[2].csat_value,b[2].csat_value)
        self.assertNotEqual(self.calc(mode='known')[0].iloc[0]['Soil Concentration Threshold [CiT] (mg/kg)'],self.calc(eff=eff,mode='known')[0].iloc[0]['Soil Concentration Threshold [CiT] (mg/kg)'])
    def test_unknown_balance(self):
        eff=app.effective_properties(self.names,self.props,dict(zip(self.names,[25.,25.])),100.)
        expected=(.25/165.85)/(.25/165.85+.25/131.4+.5/100)
        self.assertAlmostEqual(eff.iloc[0]['Mole fraction'],expected)
        self.assertLess(eff['Mole fraction'].sum(),1.)
    def test_composition_validation(self):
        with self.assertRaises(app.ValidationError):
            app.effective_properties(self.names,self.props,dict(zip(self.names,[80.,30.])),100.)
        with self.assertRaises(app.ValidationError):
            app.effective_properties([self.names[0]]*2,self.props,self.mass,100.)
    def test_miscible_requires_solubility(self):
        d=list(self.props)[-1]
        with self.assertRaises(app.ValidationError):
            app.effective_properties([d],self.props,{d:100.},100.)
        self.props[d]['solubility']=100000.
        eff=app.effective_properties([d],self.props,{d:100.},100.)
        self.assertEqual(eff.iloc[0]['Effective solubility (mg/L)'],100000.)
    def test_all_unit_checks(self):
        for index,units in [(0,'ug/kg'),(2,'fraction'),(3,'lb/ft3')]:
            raw=self.raw.copy();raw.loc[index,'units']=units
            with self.assertRaises(app.ValidationError): self.calc(raw)
    def test_non_numeric_results(self):
        for value in ['ND','<0.5','3 U','nan','inf','-1']:
            raw=self.raw.copy();raw.loc[0,'result']=value
            with self.assertRaises(app.ValidationError): self.calc(raw)
    def test_duplicates_and_missing(self):
        with self.assertRaises(app.ValidationError): self.calc(pd.concat([self.raw,self.raw.iloc[[0]]],ignore_index=True))
        with self.assertRaises(app.ValidationError): self.calc(self.raw.drop(index=2))
        with self.assertRaises(app.ValidationError): self.calc(self.raw.drop(index=0))
    def test_saturated_and_impossible_soil(self):
        for value in ['100','33.69272237196765']:
            raw=self.raw.copy();raw.loc[2,'result']=value
            with self.assertRaises(app.ValidationError): self.calc(raw)
        raw=self.raw.copy();raw.loc[3,'result']='2.8'
        with self.assertRaises(app.ValidationError): self.calc(raw)
    def test_wet_mass_conversion(self):
        raw=self.raw.copy();raw.loc[2,'result']=str(100*.2/1.2);raw.loc[6,'result']=str(100*.1/1.1)
        a=self.calc()[2];b=self.calc(raw,moisture_basis='Wet mass')[2]
        np.testing.assert_allclose(a.csat_value,b.csat_value)
    def test_cas_and_aliases(self):
        raw=self.raw.copy();raw.loc[0,'analyte']='Unusual lab name'
        self.calc(raw)
        raw.loc[0,'analyte']='PCE'
        self.calc(raw) # CAS overrides even a familiar but different name.
        self.assertEqual(app.identify('PCE',''), 'other')
        self.assertEqual(app.identify('PCE','999-99-9'), 'other')
        self.assertEqual(app.identify('% moisuter',''),'moisture')
        self.assertEqual(app.identify('cis-1,2-Dichloroethylene','156-59-2'),list(app.ANALYTE_DB)[2])
    def test_independent_sample_keys(self):
        raw=self.raw.copy();raw.loc[4:,'field_id']='SB-01-1'
        su=self.calc(raw)[2]
        self.assertEqual(len(su),2)
        self.assertNotEqual(su.iloc[0].csat_value,su.iloc[1].csat_value)
    def test_boundary_classification(self):
        name=self.names[0];props=copy.deepcopy(self.props)
        props[name]['solubility']=2.;props[name]['henry']=1.
        eff=app.effective_properties([name],props,{name:100.},100.)
        raw=self.raw.loc[self.raw.analyte.ne('TCE')].copy()
        raw.loc[raw.analyte.eq('Dry bulk density'),'result']='1'
        raw.loc[raw.analyte.eq('Percent moisture'),'result']='0'
        raw.loc[raw.analyte.eq('PCE'),'result']='1'
        summary=app.calculate(raw,[name],props,eff,particle_density=2.)[2]
        self.assertTrue(summary.csat_value.eq(1.).all())
        self.assertTrue(summary.dnapl_present.eq('Yes').all())
        raw.loc[raw.analyte.eq('PCE'),'result']='0.999999999'
        summary=app.calculate(raw,[name],props,eff,particle_density=2.)[2]
        self.assertTrue(summary.dnapl_present.eq('No – attributed to dissolve phase').all())
    def test_excel_export(self):
        th,dt,su,_,_=self.calc()
        blob=app.excel_bytes({'CSAT Components':dt[app.DETAIL],'CSAT by Sample':su})
        table=pd.read_excel(io.BytesIO(blob),sheet_name='CSAT by Sample')
        self.assertEqual(list(table),app.META+['csat_value','dnapl_present'])
        np.testing.assert_allclose(table.csat_value,su.csat_value)
        import openpyxl
        book=openpyxl.load_workbook(io.BytesIO(app.excel_bytes({'safe':pd.DataFrame({'id':['=1+1']})})))
        self.assertEqual(book['safe']['A2'].data_type,'s')
    def test_csv_and_xlsx_reads(self):
        class Upload(io.BytesIO): pass
        csv=Upload(self.raw.to_csv(index=False).encode());csv.name='test.csv'
        self.calc(app.read_input(csv))
        xls=Upload(app.excel_bytes({'data':self.raw}));xls.name='test.xlsx'
        self.calc(app.read_input(xls))
    def test_zero_mass_is_not_renormalized(self):
        eff=app.effective_properties(self.names,self.props,{n:0. for n in self.names},150.)
        self.assertTrue(eff['Mole fraction'].eq(0).all())
        self.assertTrue(self.calc(eff=eff,mode='known')[0]['Soil Concentration Threshold [CiT] (mg/kg)'].eq(0).all())

    def test_pattern_matching(self):
        for name in ['dry bulk den','DRY BULK DENSITY (g/cm3)','bulk density, dry','Dry bulk dens.','dry density']:
            self.assertEqual(app.identify(name,''),'density',name)
        for name in ['% moisture','perc moisture','percent moisture','% MOISTURE CONTENT','moisture, percent','percent moisuter']:
            self.assertEqual(app.identify(name,''),'moisture',name)
        for name in ['wet bulk density','bulk density','particle density','moisture holding capacity']:
            self.assertEqual(app.identify(name,''),'other',name)
        self.assertEqual(app.identify('dry bulk den','DBD'),'density')
        self.assertEqual(app.identify('% moisture','127-18-4'),self.names[0])
        self.assertEqual(app.identify('% moisture','999-99-9'),'other')

    def test_cas_normalization(self):
        for cas in [' 127-18-4 ', '127–18–4','127184','00127-18-4','127 - 18 - 4']:
            self.assertEqual(app.identify('Lab name',cas),self.names[0])
        self.assertEqual(app.identify('PCE','127-18-5'),'other')

    def test_unknown_has_no_composition_dependency(self):
        baseline=self.calc()[2]
        summary=app.calculate(self.raw,self.names,self.props,None,mode='unknown')[2]
        np.testing.assert_allclose(summary.csat_value,baseline.csat_value)
        self.props[self.names[0]]['mw']=None
        summary=app.calculate(self.raw,self.names,self.props,None,mode='unknown')[2]
        np.testing.assert_allclose(summary.csat_value,baseline.csat_value)

    def test_zero_composition_not_false_exceedance(self):
        eff=app.effective_properties(self.names,self.props,{n:0. for n in self.names},150.)
        table=self.calc(eff=eff,mode='known')[0]
        self.assertTrue(table['Threshold met or exceeded'].str.startswith('Not assessable').all())

    def test_known_threshold_comparison(self):
        table=self.calc(mode='known')[0]
        expected=table['Observed concentration (mg/kg)']>=table['Soil Concentration Threshold [CiT] (mg/kg)']
        self.assertEqual(list(table['Threshold met or exceeded']),['Yes' if v else 'No' for v in expected])

    def test_lab_patterns_end_to_end(self):
        raw=self.raw.copy()
        raw.loc[raw.analyte.eq('Percent moisture'),'analyte']='perc moisture'
        raw.loc[raw.analyte.eq('Dry bulk density'),'analyte']='dry bulk den'
        raw.loc[raw.cas_no.eq('127-18-4'),'analyte']='Lab PERC result'
        np.testing.assert_allclose(self.calc(raw)[2].csat_value,self.calc()[2].csat_value)

    def test_component_threshold_ratios_and_export(self):
        table=self.calc(mode='known')[0]
        expected=table['Observed concentration (mg/kg)']/table['Soil Concentration Threshold [CiT] (mg/kg)']
        np.testing.assert_allclose(table['Saturation Index'],expected)
        workbook=app.excel_bytes({'Component Thresholds':table})
        exported=pd.read_excel(io.BytesIO(workbook))
        np.testing.assert_allclose(exported['Saturation Index'],expected)
        eff=app.effective_properties(self.names,self.props,{n:0. for n in self.names},150.)
        zero=self.calc(eff=eff,mode='known')[0]
        self.assertTrue(zero['Saturation Index'].isna().all())
        exported=pd.read_excel(io.BytesIO(app.excel_bytes({'Thresholds':zero})))
        self.assertTrue(exported['Saturation Index'].isna().all())

    def test_component_ratio_boundary(self):
        name=self.names[0]; props=copy.deepcopy(self.props)
        props[name]['solubility']=2.; props[name]['henry']=1.
        eff=app.effective_properties([name],props,{name:100.},100.)
        raw=self.raw.loc[self.raw.analyte.ne('TCE')].copy()
        raw.loc[raw.analyte.eq('Dry bulk density'),'result']='1'
        raw.loc[raw.analyte.eq('Percent moisture'),'result']='0'
        for value, flag in [('0.5','No'),('1','Yes'),('1.5','Yes')]:
            raw.loc[raw.analyte.eq('PCE'),'result']=value
            th,detail,summary,_,_=app.calculate(raw,[name],props,eff,particle_density=2.,mode='known')
            np.testing.assert_allclose(th['Saturation Index'],float(value))
            self.assertTrue(th['Threshold met or exceeded'].eq(flag).all())
            self.assertTrue(detail.empty and summary.empty)

if __name__=='__main__': unittest.main()


