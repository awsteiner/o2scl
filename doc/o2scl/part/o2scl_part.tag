<?xml version='1.0' encoding='UTF-8' standalone='yes' ?>
<tagfile>
  <compound kind="file">
    <name>hdf_nucmass_io.h</name>
    <path>/Users/awsteiner/svn/osf/trunk/src/nuclei/</path>
    <filename>hdf__nucmass__io_8h</filename>
  </compound>
  <compound kind="file">
    <name>part.h</name>
    <path>/Users/awsteiner/svn/osf/trunk/src/part/</path>
    <filename>part_8h</filename>
    <class kind="class">o2scl::thermo</class>
    <class kind="class">o2scl::part</class>
  </compound>
  <compound kind="page">
    <name>part_section</name>
    <title>Particles</title>
    <filename>part_section</filename>
    <docanchor file="part_section" title="Particle data classes">part_data_sect</docanchor>
    <docanchor file="part_section" title="Classes for particle thermodynamics">part_alg_sect</docanchor>
    <docanchor file="part_section" title="Thermodynamics with derivatives">part_deriv_sect</docanchor>
    <docanchor file="part_section" title="Particle example">ex_part_sect</docanchor>
  </compound>
  <compound kind="page">
    <name>nuclei_section</name>
    <title>Nuclei and nuclear masses</title>
    <filename>nuclei_section</filename>
    <docanchor file="nuclei_section" title="Nuclear mass example">ex_nucmass_sect</docanchor>
    <docanchor file="nuclei_section" title="Nuclear mass fit example">ex_mass_fit_sect</docanchor>
  </compound>
  <compound kind="page">
    <name>partref_section</name>
    <title>Bibliography</title>
    <filename>partref_section</filename>
    <docanchor file="partref_section">Audi95</docanchor>
    <docanchor file="partref_section">Audi03</docanchor>
    <docanchor file="partref_section">Audi12</docanchor>
    <docanchor file="partref_section">Broderick00</docanchor>
    <docanchor file="partref_section">Callen</docanchor>
    <docanchor file="partref_section">Duflo95</docanchor>
    <docanchor file="partref_section">Eggleton73</docanchor>
    <docanchor file="partref_section">Goriely02</docanchor>
    <docanchor file="partref_section">Goriely07</docanchor>
    <docanchor file="partref_section">Johns96</docanchor>
    <docanchor file="partref_section">Koura00</docanchor>
    <docanchor file="partref_section">Koura05</docanchor>
    <docanchor file="partref_section">Landau</docanchor>
    <docanchor file="partref_section">Lunney03</docanchor>
    <docanchor file="partref_section">MendozaTemis10</docanchor>
    <docanchor file="partref_section">Moller95</docanchor>
    <docanchor file="partref_section">Moller97</docanchor>
    <docanchor file="partref_section">Samyn04</docanchor>
    <docanchor file="partref_section">Wang12</docanchor>
  </compound>
  <compound kind="class">
    <name>o2scl::boson</name>
    <filename>classo2scl_1_1boson.html</filename>
    <base>o2scl::part</base>
    <member kind="function">
      <type></type>
      <name>boson</name>
      <anchorfile>classo2scl_1_1boson.html</anchorfile>
      <anchor>a85ea85261f3101aaf0bc2b062d7d64d9</anchor>
      <arglist>(double mass=0.0, double dof=0.0)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>massless_calc</name>
      <anchorfile>classo2scl_1_1boson.html</anchorfile>
      <anchor>aaab47da97025c4cc96ee309bc49492b3</anchor>
      <arglist>(double temper)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual const char *</type>
      <name>type</name>
      <anchorfile>classo2scl_1_1boson.html</anchorfile>
      <anchor>ae3b965a323ff5a70e8b9ceb5590f99f4</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>part</name>
      <anchorfile>classo2scl_1_1part.html</anchorfile>
      <anchor>a44478fe390d85d29c58376ac7584c7cc</anchor>
      <arglist>(double mass=0.0, double dof=0.0)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>init</name>
      <anchorfile>classo2scl_1_1part.html</anchorfile>
      <anchor>a67198077b725549d3bc28f8d133167f5</anchor>
      <arglist>(double mass, double dof)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>anti</name>
      <anchorfile>classo2scl_1_1part.html</anchorfile>
      <anchor>a9e9986755e05325d320e885634b384a7</anchor>
      <arglist>(part &amp;ap)</arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>co</name>
      <anchorfile>classo2scl_1_1boson.html</anchorfile>
      <anchor>a328531c8684432c9bacf4b04f61eb382</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>g</name>
      <anchorfile>classo2scl_1_1part.html</anchorfile>
      <anchor>a77ad82e3af1d48c2d6687adefb12173e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>m</name>
      <anchorfile>classo2scl_1_1part.html</anchorfile>
      <anchor>abed9d7e634ed210ac2bbb0df95d5f16c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>n</name>
      <anchorfile>classo2scl_1_1part.html</anchorfile>
      <anchor>a0c750ecee7e5f5fa538029a5f87db5b5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>ed</name>
      <anchorfile>classo2scl_1_1part.html</anchorfile>
      <anchor>a916e77a026640f50a58e951ac26252e0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>pr</name>
      <anchorfile>classo2scl_1_1part.html</anchorfile>
      <anchor>ac48104230971f69a903d2596da7c14d6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>mu</name>
      <anchorfile>classo2scl_1_1part.html</anchorfile>
      <anchor>a2968bcbca72a5c640424f566b278de42</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>en</name>
      <anchorfile>classo2scl_1_1part.html</anchorfile>
      <anchor>a5709fa784aee542671055fd5bceec387</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>ms</name>
      <anchorfile>classo2scl_1_1part.html</anchorfile>
      <anchor>acc20473d1fb511e593ad0b090fec0d7c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>nu</name>
      <anchorfile>classo2scl_1_1part.html</anchorfile>
      <anchor>a7c0d7a1e3e9505d9f98c7270a16c56fb</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>inc_rest_mass</name>
      <anchorfile>classo2scl_1_1part.html</anchorfile>
      <anchor>ae2da7cbec654a535e6496daa5441bd1c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>non_interacting</name>
      <anchorfile>classo2scl_1_1part.html</anchorfile>
      <anchor>a97c8b757b23fbd10b8d3f9c813ed0754</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>o2scl::classical</name>
    <filename>classo2scl_1_1classical.html</filename>
    <member kind="function">
      <type></type>
      <name>classical</name>
      <anchorfile>classo2scl_1_1classical.html</anchorfile>
      <anchor>a8fea78b55a90e32a9b7e38877f72977d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>calc_mu</name>
      <anchorfile>classo2scl_1_1classical.html</anchorfile>
      <anchor>a87bbc697410c62f591b98b005be088ec</anchor>
      <arglist>(part &amp;p, double temper)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>calc_density</name>
      <anchorfile>classo2scl_1_1classical.html</anchorfile>
      <anchor>aeae9c49363aa1ca8b13763919ca36f8d</anchor>
      <arglist>(part &amp;p, double temper)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual const char *</type>
      <name>type</name>
      <anchorfile>classo2scl_1_1classical.html</anchorfile>
      <anchor>a5663e17909715e9cc57a579ffcb26f21</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>o2scl::eff_boson</name>
    <filename>classo2scl_1_1eff__boson.html</filename>
    <member kind="function">
      <type></type>
      <name>eff_boson</name>
      <anchorfile>classo2scl_1_1eff__boson.html</anchorfile>
      <anchor>aec7c6af74b982b7d108b7bbe8e98e2df</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>load_coefficients</name>
      <anchorfile>classo2scl_1_1eff__boson.html</anchorfile>
      <anchor>a0f288b729e8140328f4c4b225b386cc3</anchor>
      <arglist>(int ctype)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>calc_mu</name>
      <anchorfile>classo2scl_1_1eff__boson.html</anchorfile>
      <anchor>afb8b7d5a82b75b9fcdb499ee9c98282a</anchor>
      <arglist>(boson &amp;b, double temper)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>calc_density</name>
      <anchorfile>classo2scl_1_1eff__boson.html</anchorfile>
      <anchor>aa58ec1e21ada38dfdced893d2b135336</anchor>
      <arglist>(boson &amp;b, double temper)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>pair_mu</name>
      <anchorfile>classo2scl_1_1eff__boson.html</anchorfile>
      <anchor>ab2ed45ff8b5a6aaadd1932bdda372e57</anchor>
      <arglist>(boson &amp;b, double temper)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>pair_density</name>
      <anchorfile>classo2scl_1_1eff__boson.html</anchorfile>
      <anchor>a8b23b64219d8326a98d8051cf4009d8f</anchor>
      <arglist>(boson &amp;b, double temper)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_psi_root</name>
      <anchorfile>classo2scl_1_1eff__boson.html</anchorfile>
      <anchor>a15d2a7de06623ef270406d4ed2e365d3</anchor>
      <arglist>(root&lt; funct &gt; &amp;rp)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_density_mroot</name>
      <anchorfile>classo2scl_1_1eff__boson.html</anchorfile>
      <anchor>a10d6395291e89f6507f12542ff624603</anchor>
      <arglist>(mroot&lt; mm_funct&lt;&gt;, boost::numeric::ublas::vector&lt; double &gt;, jac_funct&lt;&gt; &gt; &amp;rp)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_meth2_root</name>
      <anchorfile>classo2scl_1_1eff__boson.html</anchorfile>
      <anchor>a9530331468b753d181d5e7433d1ca033</anchor>
      <arglist>(root&lt; funct &gt; &amp;rp)</arglist>
    </member>
    <member kind="variable">
      <type>mroot_hybrids&lt; mm_funct&lt;&gt;, boost::numeric::ublas::vector&lt; double &gt;, boost::numeric::ublas::matrix&lt; double &gt;, jac_funct&lt;&gt; &gt;</type>
      <name>def_density_mroot</name>
      <anchorfile>classo2scl_1_1eff__boson.html</anchorfile>
      <anchor>a590f54bb7e700963e007539756064dbb</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>root_cern&lt; funct &gt;</type>
      <name>def_psi_root</name>
      <anchorfile>classo2scl_1_1eff__boson.html</anchorfile>
      <anchor>a064b35931453742dbd53d2a9540d7a43</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>root_cern&lt; funct &gt;</type>
      <name>def_meth2_root</name>
      <anchorfile>classo2scl_1_1eff__boson.html</anchorfile>
      <anchor>a13e06b8ff07be34ead297366b3cbf0f5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const int</type>
      <name>cf_boselat3</name>
      <anchorfile>classo2scl_1_1eff__boson.html</anchorfile>
      <anchor>a3be52c7bc1cc55550258622ab33ae0a7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const int</type>
      <name>cf_bosejel21</name>
      <anchorfile>classo2scl_1_1eff__boson.html</anchorfile>
      <anchor>aab5433d69d159f0de926f6bd559bd4f8</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const int</type>
      <name>cf_bosejel22</name>
      <anchorfile>classo2scl_1_1eff__boson.html</anchorfile>
      <anchor>a5311197232a36c1fd99942ef9483798c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const int</type>
      <name>cf_bosejel34</name>
      <anchorfile>classo2scl_1_1eff__boson.html</anchorfile>
      <anchor>a6b8ff3ffcf4435f4fade79c5af8f63b6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const int</type>
      <name>cf_bosejel34cons</name>
      <anchorfile>classo2scl_1_1eff__boson.html</anchorfile>
      <anchor>aaf7dfc410c80c75733cb42e1622ba32a</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" protection="protected">
      <type>double</type>
      <name>solve_fun</name>
      <anchorfile>classo2scl_1_1eff__boson.html</anchorfile>
      <anchor>aff436e99b721db79661817c2a166505a</anchor>
      <arglist>(double x, double &amp;psi)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>int</type>
      <name>density_fun</name>
      <anchorfile>classo2scl_1_1eff__boson.html</anchorfile>
      <anchor>a126078530d06467ee2088c4c0c447e38</anchor>
      <arglist>(size_t nv, const ubvector &amp;x, ubvector &amp;y)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>int</type>
      <name>pair_density_fun</name>
      <anchorfile>classo2scl_1_1eff__boson.html</anchorfile>
      <anchor>a62b54cd40c52c02c9772dbb9c812827c</anchor>
      <arglist>(size_t nv, const ubvector &amp;x, ubvector &amp;y)</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>ubmatrix</type>
      <name>Pmnb</name>
      <anchorfile>classo2scl_1_1eff__boson.html</anchorfile>
      <anchor>a5637a12137015826abf89a8844883fa4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>int</type>
      <name>sizem</name>
      <anchorfile>classo2scl_1_1eff__boson.html</anchorfile>
      <anchor>af086f5e4c344cea8ab4d054901da28bc</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>int</type>
      <name>sizen</name>
      <anchorfile>classo2scl_1_1eff__boson.html</anchorfile>
      <anchor>a3892e0cb955c2c715426e1169c1d8114</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>double</type>
      <name>parma</name>
      <anchorfile>classo2scl_1_1eff__boson.html</anchorfile>
      <anchor>a34229294e57ae67581479f09faea07f4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>double</type>
      <name>fix_density</name>
      <anchorfile>classo2scl_1_1eff__boson.html</anchorfile>
      <anchor>a7b61f16c0d246f685e752b3f6bfcbd2b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>boson *</type>
      <name>bp</name>
      <anchorfile>classo2scl_1_1eff__boson.html</anchorfile>
      <anchor>a5237db2b281be0607fefbb4815e3a20e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>double</type>
      <name>T</name>
      <anchorfile>classo2scl_1_1eff__boson.html</anchorfile>
      <anchor>ac1f10471170ac0944513c6670e46d7a9</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>mroot&lt; mm_funct&lt;&gt;, boost::numeric::ublas::vector&lt; double &gt;, jac_funct&lt;&gt; &gt; *</type>
      <name>density_mroot</name>
      <anchorfile>classo2scl_1_1eff__boson.html</anchorfile>
      <anchor>a9a3cae9a1e0b653e9e19c32a3e8baef6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>root&lt; funct &gt; *</type>
      <name>psi_root</name>
      <anchorfile>classo2scl_1_1eff__boson.html</anchorfile>
      <anchor>a98ebe5a93451f5d4f5bf1d1d4c02d829</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>root&lt; funct &gt; *</type>
      <name>meth2_root</name>
      <anchorfile>classo2scl_1_1eff__boson.html</anchorfile>
      <anchor>aba1301dcb5a7d148b4efd1e3bc081b80</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>o2scl::eff_fermion</name>
    <filename>classo2scl_1_1eff__fermion.html</filename>
    <base>o2scl::fermion_eval_thermo</base>
    <class kind="class">o2scl::eff_fermion::density_fun</class>
    <class kind="class">o2scl::eff_fermion::pair_density_fun</class>
    <member kind="function">
      <type></type>
      <name>eff_fermion</name>
      <anchorfile>classo2scl_1_1eff__fermion.html</anchorfile>
      <anchor>a08cd91f4943f8037abbe16c51d91abae</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>calc_mu</name>
      <anchorfile>classo2scl_1_1eff__fermion.html</anchorfile>
      <anchor>a5a6710e214d9253fdbf4ee5f45f9bc6a</anchor>
      <arglist>(fermion &amp;f, double temper)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>calc_density</name>
      <anchorfile>classo2scl_1_1eff__fermion.html</anchorfile>
      <anchor>a1c0a663484cf1da37a04f5e34cd25e4c</anchor>
      <arglist>(fermion &amp;f, double temper)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>pair_mu</name>
      <anchorfile>classo2scl_1_1eff__fermion.html</anchorfile>
      <anchor>ab4e4307826f6d971f54457af48343b9c</anchor>
      <arglist>(fermion &amp;f, double temper)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>pair_density</name>
      <anchorfile>classo2scl_1_1eff__fermion.html</anchorfile>
      <anchor>a2b1cdf061ffed49fd6b3d0d4d11e5e6c</anchor>
      <arglist>(fermion &amp;f, double temper)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>set_psi_root</name>
      <anchorfile>classo2scl_1_1eff__fermion.html</anchorfile>
      <anchor>ab5cdb3278a3b88fc531b4e71a893a25e</anchor>
      <arglist>(root&lt; funct &gt; &amp;rp)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>set_density_root</name>
      <anchorfile>classo2scl_1_1eff__fermion.html</anchorfile>
      <anchor>a93efb087678af511d2f8b0c61985ef90</anchor>
      <arglist>(root&lt; funct &gt; &amp;rp)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual const char *</type>
      <name>type</name>
      <anchorfile>classo2scl_1_1eff__fermion.html</anchorfile>
      <anchor>a2d6867fd8531440902c5cc80afea8ebe</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>calc_mu_ndeg</name>
      <anchorfile>classo2scl_1_1fermion__eval__thermo.html</anchorfile>
      <anchor>a42bdbb561725a41710d6a55573e669ed</anchor>
      <arglist>(fermion &amp;f, double temper, double prec=1.0e-18)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>calc_mu_deg</name>
      <anchorfile>classo2scl_1_1fermion__eval__thermo.html</anchorfile>
      <anchor>af1110f365785dec1d317f85731677c67</anchor>
      <arglist>(fermion &amp;f, double temper, double prec=1.0e-18)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_massless_root</name>
      <anchorfile>classo2scl_1_1fermion__eval__thermo.html</anchorfile>
      <anchor>ae9ecd640a3d6bda06ca628d39a3934fc</anchor>
      <arglist>(root&lt; funct &gt; &amp;rp)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>calibrate</name>
      <anchorfile>classo2scl_1_1fermion__eval__thermo.html</anchorfile>
      <anchor>a6c4d655ea6b0aa83556c0e2173ece649</anchor>
      <arglist>(fermion &amp;f, int verbose=0, std::string fname=&quot;&quot;)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>massless_calc_mu</name>
      <anchorfile>classo2scl_1_1fermion__eval__thermo.html</anchorfile>
      <anchor>a3c101369352f14df09d9acf48fd41a9b</anchor>
      <arglist>(fermion &amp;f, double temper)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>massless_calc_density</name>
      <anchorfile>classo2scl_1_1fermion__eval__thermo.html</anchorfile>
      <anchor>a17a03a6bd387f1ba15d1985d6a63bdc9</anchor>
      <arglist>(fermion &amp;f, double temper)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>massless_pair_mu</name>
      <anchorfile>classo2scl_1_1fermion__eval__thermo.html</anchorfile>
      <anchor>a94837ae3fc57184ae938d868e1eb85c0</anchor>
      <arglist>(fermion &amp;f, double temper)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>massless_pair_density</name>
      <anchorfile>classo2scl_1_1fermion__eval__thermo.html</anchorfile>
      <anchor>a349fb84683acb4eeca137ff3ee1b0ebc</anchor>
      <arglist>(fermion &amp;f, double temper)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>kf_from_density</name>
      <anchorfile>classo2scl_1_1fermion__zerot.html</anchorfile>
      <anchor>a5f333606949ef1cf6318fefce83ec87e</anchor>
      <arglist>(fermion &amp;f)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>energy_density_zerot</name>
      <anchorfile>classo2scl_1_1fermion__zerot.html</anchorfile>
      <anchor>aad7ae2e73c436bb2f0c302e073c5e376</anchor>
      <arglist>(fermion &amp;f)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>pressure_zerot</name>
      <anchorfile>classo2scl_1_1fermion__zerot.html</anchorfile>
      <anchor>a127f543719bae19673248f4efc9ef806</anchor>
      <arglist>(fermion &amp;f)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>calc_mu_zerot</name>
      <anchorfile>classo2scl_1_1fermion__zerot.html</anchorfile>
      <anchor>ac2eb379599c918bb04abaddc8ead587f</anchor>
      <arglist>(fermion &amp;f)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>calc_density_zerot</name>
      <anchorfile>classo2scl_1_1fermion__zerot.html</anchorfile>
      <anchor>a11324827b7ae7a32a633851e8908d9bb</anchor>
      <arglist>(fermion &amp;f)</arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>tlimit</name>
      <anchorfile>classo2scl_1_1eff__fermion.html</anchorfile>
      <anchor>a7775655adb09365b806cf00ab77aadd7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>root_cern&lt; funct &gt;</type>
      <name>def_psi_root</name>
      <anchorfile>classo2scl_1_1eff__fermion.html</anchorfile>
      <anchor>aa6011739ae61ce67260abab1652c1dff</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>root_cern&lt; funct &gt;</type>
      <name>def_density_root</name>
      <anchorfile>classo2scl_1_1eff__fermion.html</anchorfile>
      <anchor>a8bbaa3fdbaea11f847c1b3939f08b388</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>min_psi</name>
      <anchorfile>classo2scl_1_1eff__fermion.html</anchorfile>
      <anchor>a9acf4ae1344255e382c1f45a504dcf30</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>root_cern&lt; funct &gt;</type>
      <name>def_massless_root</name>
      <anchorfile>classo2scl_1_1fermion__eval__thermo.html</anchorfile>
      <anchor>a585ab38f871f62a53848cab2f0344295</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" protection="protected">
      <type>double</type>
      <name>solve_fun</name>
      <anchorfile>classo2scl_1_1eff__fermion.html</anchorfile>
      <anchor>addb294b7fd4fbe5658dfa3777acdc8c3</anchor>
      <arglist>(double x, double &amp;psi)</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>ubmatrix</type>
      <name>Pmnf</name>
      <anchorfile>classo2scl_1_1eff__fermion.html</anchorfile>
      <anchor>aa28a898d9c2dbdf3a5faa492cffabb98</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>double</type>
      <name>parma</name>
      <anchorfile>classo2scl_1_1eff__fermion.html</anchorfile>
      <anchor>a3c5987c49fbb1e4a43a53c778a94fb97</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>int</type>
      <name>sizem</name>
      <anchorfile>classo2scl_1_1eff__fermion.html</anchorfile>
      <anchor>a57544338f6f02c8c9ea4f18ee5900697</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>int</type>
      <name>sizen</name>
      <anchorfile>classo2scl_1_1eff__fermion.html</anchorfile>
      <anchor>aca226229e1f4ad496c18934c50852eba</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>root&lt; funct &gt; *</type>
      <name>psi_root</name>
      <anchorfile>classo2scl_1_1eff__fermion.html</anchorfile>
      <anchor>a9974f213a86c6c60195dc79a37b7ae42</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>root&lt; funct &gt; *</type>
      <name>density_root</name>
      <anchorfile>classo2scl_1_1eff__fermion.html</anchorfile>
      <anchor>ac32af4f86163389e2a276840c35bb1c6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>root&lt; funct &gt; *</type>
      <name>massless_root</name>
      <anchorfile>classo2scl_1_1fermion__eval__thermo.html</anchorfile>
      <anchor>a149161e3745807926b08b1a9905e1f03</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const int</type>
      <name>cf_fermilat3</name>
      <anchorfile>classo2scl_1_1eff__fermion.html</anchorfile>
      <anchor>a750e6a26665515824a4cc4bc02accb05</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const int</type>
      <name>cf_fermijel2</name>
      <anchorfile>classo2scl_1_1eff__fermion.html</anchorfile>
      <anchor>ac0d25713b6316d0a35095dc8c5a3db74</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const int</type>
      <name>cf_fermijel3</name>
      <anchorfile>classo2scl_1_1eff__fermion.html</anchorfile>
      <anchor>ab2eb508cfe45008d9e2f39aeea723d76</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const int</type>
      <name>cf_fermijel3cons</name>
      <anchorfile>classo2scl_1_1eff__fermion.html</anchorfile>
      <anchor>a946538db5f597eca297a84515b543469</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>load_coefficients</name>
      <anchorfile>classo2scl_1_1eff__fermion.html</anchorfile>
      <anchor>a8cc886446c330d92a1bdbb08e8a801c5</anchor>
      <arglist>(int ctype)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>o2scl::eff_fermion::density_fun</name>
    <filename>classo2scl_1_1eff__fermion_1_1density__fun.html</filename>
    <base>o2scl::funct</base>
    <member kind="function">
      <type>double</type>
      <name>operator()</name>
      <anchorfile>classo2scl_1_1eff__fermion_1_1density__fun.html</anchorfile>
      <anchor>a4ba9c59e81564cdc6fe669ca7a3f06bc</anchor>
      <arglist>(double x) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>o2scl::eff_fermion::pair_density_fun</name>
    <filename>classo2scl_1_1eff__fermion_1_1pair__density__fun.html</filename>
    <base>o2scl::funct</base>
    <member kind="function">
      <type>double</type>
      <name>operator()</name>
      <anchorfile>classo2scl_1_1eff__fermion_1_1pair__density__fun.html</anchorfile>
      <anchor>a18e005c7185bba5c853458fa1a560027</anchor>
      <arglist>(double x) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>o2scl::fermion</name>
    <filename>classo2scl_1_1fermion.html</filename>
    <base>o2scl::part</base>
    <member kind="function">
      <type></type>
      <name>fermion</name>
      <anchorfile>classo2scl_1_1fermion.html</anchorfile>
      <anchor>ae5b585a4d06fd84fefd985e256f64b66</anchor>
      <arglist>(double mass=0, double dof=0)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual const char *</type>
      <name>type</name>
      <anchorfile>classo2scl_1_1fermion.html</anchorfile>
      <anchor>a17425c94820ee667a2cca3d5f47c80fe</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>kf</name>
      <anchorfile>classo2scl_1_1fermion.html</anchorfile>
      <anchor>af508160a9c3b06cb7014f977b5ad3e70</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>del</name>
      <anchorfile>classo2scl_1_1fermion.html</anchorfile>
      <anchor>aa8e723969acd716dab597598a090bad5</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>o2scl::fermion_zerot</name>
    <filename>classo2scl_1_1fermion__zerot.html</filename>
  </compound>
  <compound kind="class">
    <name>o2scl::fermion_eval_thermo</name>
    <filename>classo2scl_1_1fermion__eval__thermo.html</filename>
    <base>o2scl::fermion_zerot</base>
    <class kind="class">o2scl::fermion_eval_thermo::massless_fun</class>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>calc_mu</name>
      <anchorfile>classo2scl_1_1fermion__eval__thermo.html</anchorfile>
      <anchor>a6016bf191e0ec968d4c7979f1469127a</anchor>
      <arglist>(fermion &amp;f, double temper)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>calc_density</name>
      <anchorfile>classo2scl_1_1fermion__eval__thermo.html</anchorfile>
      <anchor>ad958390199ace5bbda7b9a8738e17ff7</anchor>
      <arglist>(fermion &amp;f, double temper)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>pair_mu</name>
      <anchorfile>classo2scl_1_1fermion__eval__thermo.html</anchorfile>
      <anchor>a8a5a2d7bb77ce182262e46f11fd3b87e</anchor>
      <arglist>(fermion &amp;f, double temper)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>pair_density</name>
      <anchorfile>classo2scl_1_1fermion__eval__thermo.html</anchorfile>
      <anchor>a074fab18aeb08e6d0afd3a20f69d6829</anchor>
      <arglist>(fermion &amp;f, double temper)=0</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual const char *</type>
      <name>type</name>
      <anchorfile>classo2scl_1_1fermion__eval__thermo.html</anchorfile>
      <anchor>a6418b55eaa456511d94d7fa2a40714d5</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>o2scl::fermion_eval_thermo::massless_fun</name>
    <filename>classo2scl_1_1fermion__eval__thermo_1_1massless__fun.html</filename>
    <base>o2scl::funct</base>
    <member kind="function">
      <type>double</type>
      <name>operator()</name>
      <anchorfile>classo2scl_1_1fermion__eval__thermo_1_1massless__fun.html</anchorfile>
      <anchor>a7e0b06463e14cab7d7d743cf5ee4dc1c</anchor>
      <arglist>(double x) const </arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>fermion &amp;</type>
      <name>f_</name>
      <anchorfile>classo2scl_1_1fermion__eval__thermo_1_1massless__fun.html</anchorfile>
      <anchor>a28352befd29ec30cdf97404a2e74711c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>double</type>
      <name>temper_</name>
      <anchorfile>classo2scl_1_1fermion__eval__thermo_1_1massless__fun.html</anchorfile>
      <anchor>a8c1bc745b0e3d1a02213eec039682b4d</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>o2scl::mag_fermion_zerot</name>
    <filename>classo2scl_1_1mag__fermion__zerot.html</filename>
    <base>o2scl::fermion_zerot</base>
    <member kind="function">
      <type></type>
      <name>mag_fermion_zerot</name>
      <anchorfile>classo2scl_1_1mag__fermion__zerot.html</anchorfile>
      <anchor>a891511b5df9488978a7a02506ac0f160</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>calc_mu_zerot_mag</name>
      <anchorfile>classo2scl_1_1mag__fermion__zerot.html</anchorfile>
      <anchor>acc92160fefa3235e16785428af185eb7</anchor>
      <arglist>(fermion &amp;f, double qB, double kappa=0.0)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>calc_density_zerot_mag</name>
      <anchorfile>classo2scl_1_1mag__fermion__zerot.html</anchorfile>
      <anchor>a23b2024cd5062dca12b1c3e3900e8c55</anchor>
      <arglist>(fermion &amp;f, double qB, double kappa=0.0)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>set_density_root</name>
      <anchorfile>classo2scl_1_1mag__fermion__zerot.html</anchorfile>
      <anchor>ac231ab663e71b890d619151682c105eb</anchor>
      <arglist>(mroot&lt; mm_funct&lt;&gt;, boost::numeric::ublas::vector&lt; double &gt;, jac_funct&lt;&gt; &gt; &amp;rp)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual const char *</type>
      <name>type</name>
      <anchorfile>classo2scl_1_1mag__fermion__zerot.html</anchorfile>
      <anchor>a9e684164dd609a7fc3e8f37b9a6c1494</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>nmax_up</name>
      <anchorfile>classo2scl_1_1mag__fermion__zerot.html</anchorfile>
      <anchor>a43afd80d6ed530a7cf41b322495ca869</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>nmax_dn</name>
      <anchorfile>classo2scl_1_1mag__fermion__zerot.html</anchorfile>
      <anchor>a913e94672131f6a31f0176bd5d334d26</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>sum_limit</name>
      <anchorfile>classo2scl_1_1mag__fermion__zerot.html</anchorfile>
      <anchor>a13318f4fa0a7c1e64ee1ddfd83176924</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>mroot_hybrids&lt; mm_funct&lt;&gt;, boost::numeric::ublas::vector&lt; double &gt;, boost::numeric::ublas::matrix&lt; double &gt;, jac_funct&lt;&gt; &gt;</type>
      <name>def_density_root</name>
      <anchorfile>classo2scl_1_1mag__fermion__zerot.html</anchorfile>
      <anchor>ab8efc684591d0499f8af01417b63d1cf</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" protection="protected">
      <type>int</type>
      <name>solve_fun</name>
      <anchorfile>classo2scl_1_1mag__fermion__zerot.html</anchorfile>
      <anchor>aae2ffc7d9521bd2cf38a4fb0f374fbda</anchor>
      <arglist>(size_t nv, const ubvector &amp;x, ubvector &amp;y)</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>mroot&lt; mm_funct&lt;&gt;, boost::numeric::ublas::vector&lt; double &gt;, jac_funct&lt;&gt; &gt; *</type>
      <name>density_root</name>
      <anchorfile>classo2scl_1_1mag__fermion__zerot.html</anchorfile>
      <anchor>a024d46652a869472192d5a0240f738e9</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>double</type>
      <name>qBt</name>
      <anchorfile>classo2scl_1_1mag__fermion__zerot.html</anchorfile>
      <anchor>ab205ad90d2afb40c9b002ce85376cf98</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>double</type>
      <name>kt</name>
      <anchorfile>classo2scl_1_1mag__fermion__zerot.html</anchorfile>
      <anchor>ab6b2817c7c0bd14ec2d1516c0e66caad</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>double</type>
      <name>dent</name>
      <anchorfile>classo2scl_1_1mag__fermion__zerot.html</anchorfile>
      <anchor>a3029dbdb6112777e5e687ef28b4fa2a9</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>fermion *</type>
      <name>fp</name>
      <anchorfile>classo2scl_1_1mag__fermion__zerot.html</anchorfile>
      <anchor>a48d3533fd919fbde51a3de15c19edc3f</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>o2scl::nonrel_fermion</name>
    <filename>classo2scl_1_1nonrel__fermion.html</filename>
    <base>o2scl::fermion_eval_thermo</base>
    <member kind="function">
      <type></type>
      <name>nonrel_fermion</name>
      <anchorfile>classo2scl_1_1nonrel__fermion.html</anchorfile>
      <anchor>a560c80b93f12cb11af315046464326b4</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>calc_mu_zerot</name>
      <anchorfile>classo2scl_1_1nonrel__fermion.html</anchorfile>
      <anchor>ae6b6e8d32e45490ffad2be5075d90ecb</anchor>
      <arglist>(fermion &amp;f)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>calc_density_zerot</name>
      <anchorfile>classo2scl_1_1nonrel__fermion.html</anchorfile>
      <anchor>a1f50a5cd1e76ad56807a948aad9e0a1e</anchor>
      <arglist>(fermion &amp;f)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>calc_mu</name>
      <anchorfile>classo2scl_1_1nonrel__fermion.html</anchorfile>
      <anchor>a2c142746a9d40b005239a0073f72eeff</anchor>
      <arglist>(fermion &amp;f, double temper)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>calc_density</name>
      <anchorfile>classo2scl_1_1nonrel__fermion.html</anchorfile>
      <anchor>abbb0e7cfb69a56af54f20afe922df4b2</anchor>
      <arglist>(fermion &amp;f, double temper)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>pair_mu</name>
      <anchorfile>classo2scl_1_1nonrel__fermion.html</anchorfile>
      <anchor>ac66ace49f48c277431ce3866328ea883</anchor>
      <arglist>(fermion &amp;f, double temper)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>pair_density</name>
      <anchorfile>classo2scl_1_1nonrel__fermion.html</anchorfile>
      <anchor>a92a46d6c5d072b55a8b5df950b1ac056</anchor>
      <arglist>(fermion &amp;f, double temper)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>nu_from_n</name>
      <anchorfile>classo2scl_1_1nonrel__fermion.html</anchorfile>
      <anchor>adc41da8e5a5c913a8410491a20231ca4</anchor>
      <arglist>(fermion &amp;f, double temper)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>set_density_root</name>
      <anchorfile>classo2scl_1_1nonrel__fermion.html</anchorfile>
      <anchor>a67013e3a91c31bef0d11f3822a27ac34</anchor>
      <arglist>(root&lt; funct &gt; &amp;rp)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual const char *</type>
      <name>type</name>
      <anchorfile>classo2scl_1_1nonrel__fermion.html</anchorfile>
      <anchor>a7c7924184baa2e45cf4dcf016d78f554</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable">
      <type>root_cern&lt; funct &gt;</type>
      <name>def_density_root</name>
      <anchorfile>classo2scl_1_1nonrel__fermion.html</anchorfile>
      <anchor>a34cff9eac8dace35a7608ddb52bb0460</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" protection="protected">
      <type>double</type>
      <name>solve_fun</name>
      <anchorfile>classo2scl_1_1nonrel__fermion.html</anchorfile>
      <anchor>a643c82e5029de0268c42bd41b4185b19</anchor>
      <arglist>(double x)</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>root&lt; funct &gt; *</type>
      <name>density_root</name>
      <anchorfile>classo2scl_1_1nonrel__fermion.html</anchorfile>
      <anchor>aa45c4247d402d9007bd4129ad2ad8f51</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>o2scl::thermo</name>
    <filename>classo2scl_1_1thermo.html</filename>
    <member kind="function">
      <type>const char *</type>
      <name>type</name>
      <anchorfile>classo2scl_1_1thermo.html</anchorfile>
      <anchor>a62b63ff8533c549ba91d4cb565d02dd8</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>pr</name>
      <anchorfile>classo2scl_1_1thermo.html</anchorfile>
      <anchor>aaefa30b0eb7d0987b399db4737e9b9e1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>ed</name>
      <anchorfile>classo2scl_1_1thermo.html</anchorfile>
      <anchor>a4c3156be062fc1df3abae766b9a22a01</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>en</name>
      <anchorfile>classo2scl_1_1thermo.html</anchorfile>
      <anchor>ac4e695c90719b3347d5f5f5c0077cda2</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>o2scl::part</name>
    <filename>classo2scl_1_1part.html</filename>
    <member kind="function" virtualness="virtual">
      <type>virtual const char *</type>
      <name>type</name>
      <anchorfile>classo2scl_1_1part.html</anchorfile>
      <anchor>a609db0ef551e11d098816e41a8e6eb00</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>o2scl::part_deriv</name>
    <filename>classo2scl_1_1part__deriv.html</filename>
    <base>o2scl::part</base>
    <member kind="function">
      <type></type>
      <name>part_deriv</name>
      <anchorfile>classo2scl_1_1part__deriv.html</anchorfile>
      <anchor>a5af2cb28303c23245c683b516efe26c9</anchor>
      <arglist>(double mass=0.0, double dof=0.0)</arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>dndmu</name>
      <anchorfile>classo2scl_1_1part__deriv.html</anchorfile>
      <anchor>afb71dc0012ba41645817b9d4fde3823a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>dndT</name>
      <anchorfile>classo2scl_1_1part__deriv.html</anchorfile>
      <anchor>af692ba6c8a28981163d671766d92f7fc</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>dsdT</name>
      <anchorfile>classo2scl_1_1part__deriv.html</anchorfile>
      <anchor>a0bc07a75aa4b65dc32eca32895012815</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>dndm</name>
      <anchorfile>classo2scl_1_1part__deriv.html</anchorfile>
      <anchor>a6122dd70e7be2a9b4c87731475e3bfb7</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>o2scl::fermion_deriv</name>
    <filename>classo2scl_1_1fermion__deriv.html</filename>
    <base>o2scl::part_deriv</base>
    <member kind="function">
      <type></type>
      <name>fermion_deriv</name>
      <anchorfile>classo2scl_1_1fermion__deriv.html</anchorfile>
      <anchor>a9668b9dbe6d1950b0025c7febf7c141e</anchor>
      <arglist>(double mass=0.0, double dof=0.0)</arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>kf</name>
      <anchorfile>classo2scl_1_1fermion__deriv.html</anchorfile>
      <anchor>ae15b08b3e84437d2be160a7bbd15ac3b</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>o2scl::fermion_deriv_thermo</name>
    <filename>classo2scl_1_1fermion__deriv__thermo.html</filename>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>calc_mu</name>
      <anchorfile>classo2scl_1_1fermion__deriv__thermo.html</anchorfile>
      <anchor>a51c00c7c4782b05c15d4254056c9fbc4</anchor>
      <arglist>(fermion_deriv &amp;f, double temper)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>calc_density</name>
      <anchorfile>classo2scl_1_1fermion__deriv__thermo.html</anchorfile>
      <anchor>aee885cf624cc5d4f55167df9229dde33</anchor>
      <arglist>(fermion_deriv &amp;f, double temper)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>pair_mu</name>
      <anchorfile>classo2scl_1_1fermion__deriv__thermo.html</anchorfile>
      <anchor>a59e48bdf84d692057eaeb26515a3c607</anchor>
      <arglist>(fermion_deriv &amp;f, double temper)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>pair_density</name>
      <anchorfile>classo2scl_1_1fermion__deriv__thermo.html</anchorfile>
      <anchor>a3033fe1fd608e89f9af5e0ecbac52edd</anchor>
      <arglist>(fermion_deriv &amp;f, double temper)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>nu_from_n</name>
      <anchorfile>classo2scl_1_1fermion__deriv__thermo.html</anchorfile>
      <anchor>ac5a761b91b352b55b171280595006114</anchor>
      <arglist>(fermion_deriv &amp;f, double temper)=0</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>o2scl::quark</name>
    <filename>classo2scl_1_1quark.html</filename>
    <base>o2scl::fermion</base>
    <member kind="function">
      <type></type>
      <name>quark</name>
      <anchorfile>classo2scl_1_1quark.html</anchorfile>
      <anchor>a9d7b06b3df1d5507526604b966ff8329</anchor>
      <arglist>(double m=0.0, double g=0.0)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual const char *</type>
      <name>type</name>
      <anchorfile>classo2scl_1_1quark.html</anchorfile>
      <anchor>a59a852ec4d59ec95bbd794068a12f2d5</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>B</name>
      <anchorfile>classo2scl_1_1quark.html</anchorfile>
      <anchor>a11bf4ddd40b05be4200b0a3bb6f3a1fc</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>qq</name>
      <anchorfile>classo2scl_1_1quark.html</anchorfile>
      <anchor>a4200983215bcf52bfc0b1fc2153c131d</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>o2scl::rel_boson</name>
    <filename>classo2scl_1_1rel__boson.html</filename>
    <member kind="function">
      <type></type>
      <name>rel_boson</name>
      <anchorfile>classo2scl_1_1rel__boson.html</anchorfile>
      <anchor>aa993f0ed59662e3811ad04975421a2dc</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>calc_mu</name>
      <anchorfile>classo2scl_1_1rel__boson.html</anchorfile>
      <anchor>a2c06a8299f29f3c694b5cc85f424ae61</anchor>
      <arglist>(boson &amp;b, double temper)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>calc_density</name>
      <anchorfile>classo2scl_1_1rel__boson.html</anchorfile>
      <anchor>a35884d0f58afbf564c901ba91503c178</anchor>
      <arglist>(boson &amp;b, double temper)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>pair_mu</name>
      <anchorfile>classo2scl_1_1rel__boson.html</anchorfile>
      <anchor>ac258331292e6938ac0da526c0ac4f749</anchor>
      <arglist>(boson &amp;b, double temper)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>pair_density</name>
      <anchorfile>classo2scl_1_1rel__boson.html</anchorfile>
      <anchor>a981c11fd296e5acbb7134f4d857e4a01</anchor>
      <arglist>(boson &amp;b, double temper)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>nu_from_n</name>
      <anchorfile>classo2scl_1_1rel__boson.html</anchorfile>
      <anchor>ade1c7aac7af641af1ab1745e94134fab</anchor>
      <arglist>(boson &amp;b, double temper)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_inte</name>
      <anchorfile>classo2scl_1_1rel__boson.html</anchorfile>
      <anchor>a954f53b031f4ea8ab2284968adae575e</anchor>
      <arglist>(inte&lt; funct &gt; &amp;l_nit, inte&lt; funct &gt; &amp;l_dit)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_density_root</name>
      <anchorfile>classo2scl_1_1rel__boson.html</anchorfile>
      <anchor>aade822d1ee7207c23b61ac8e102f537d</anchor>
      <arglist>(root&lt; funct &gt; &amp;rp)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual const char *</type>
      <name>type</name>
      <anchorfile>classo2scl_1_1rel__boson.html</anchorfile>
      <anchor>a727a2782d967744f7ecb1d75b8856ecd</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>mroot_err</name>
      <anchorfile>classo2scl_1_1rel__boson.html</anchorfile>
      <anchor>af8fc5f6ddb62dbfc4507f081deb2fae7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>inte_err</name>
      <anchorfile>classo2scl_1_1rel__boson.html</anchorfile>
      <anchor>a7cb0aae1cc0ecba447e8e12c16501611</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>root_cern&lt; funct &gt;</type>
      <name>def_density_root</name>
      <anchorfile>classo2scl_1_1rel__boson.html</anchorfile>
      <anchor>a8438a41532b90a6f33ceb43d93bdb364</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>inte_qagiu_gsl&lt; funct &gt;</type>
      <name>def_nit</name>
      <anchorfile>classo2scl_1_1rel__boson.html</anchorfile>
      <anchor>a56f9de66c2d1aba8904959c21107fb6f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>inte_qag_gsl&lt; funct &gt;</type>
      <name>def_dit</name>
      <anchorfile>classo2scl_1_1rel__boson.html</anchorfile>
      <anchor>a80b4b0a81e177df7b562b7769279a6f0</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" protection="protected">
      <type>double</type>
      <name>density_fun</name>
      <anchorfile>classo2scl_1_1rel__boson.html</anchorfile>
      <anchor>ac5bf6d73cca371ac4903fd4995ce4842</anchor>
      <arglist>(double u)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>double</type>
      <name>energy_fun</name>
      <anchorfile>classo2scl_1_1rel__boson.html</anchorfile>
      <anchor>a14ead352870a9f1c66917fac9019346f</anchor>
      <arglist>(double u)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>double</type>
      <name>entropy_fun</name>
      <anchorfile>classo2scl_1_1rel__boson.html</anchorfile>
      <anchor>a8da7930016c9874965d660d0797d75d0</anchor>
      <arglist>(double u)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>double</type>
      <name>deg_density_fun</name>
      <anchorfile>classo2scl_1_1rel__boson.html</anchorfile>
      <anchor>a15bc0c0462ff0aef59198e6ae5eb6901</anchor>
      <arglist>(double u)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>double</type>
      <name>deg_energy_fun</name>
      <anchorfile>classo2scl_1_1rel__boson.html</anchorfile>
      <anchor>a916d4d669e5980fb14635abb678390e4</anchor>
      <arglist>(double u)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>double</type>
      <name>deg_entropy_fun</name>
      <anchorfile>classo2scl_1_1rel__boson.html</anchorfile>
      <anchor>a4800247b83e3b2ce5251eaada75860fc</anchor>
      <arglist>(double u)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>double</type>
      <name>solve_fun</name>
      <anchorfile>classo2scl_1_1rel__boson.html</anchorfile>
      <anchor>acc4f27f0f0c9fd321eda8632c46583f8</anchor>
      <arglist>(double x)</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>inte&lt; funct &gt; *</type>
      <name>nit</name>
      <anchorfile>classo2scl_1_1rel__boson.html</anchorfile>
      <anchor>ae4355e27ec2270a7422b7bbb2e3fa326</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>inte&lt; funct &gt; *</type>
      <name>dit</name>
      <anchorfile>classo2scl_1_1rel__boson.html</anchorfile>
      <anchor>ade0783572aa17a810dc5d01674543f51</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>root&lt; funct &gt; *</type>
      <name>density_root</name>
      <anchorfile>classo2scl_1_1rel__boson.html</anchorfile>
      <anchor>a5592215b38510467c0732cbbdf3ef312</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>boson *</type>
      <name>bp</name>
      <anchorfile>classo2scl_1_1rel__boson.html</anchorfile>
      <anchor>a916e8ef6c84a6c67a64247267bfcb365</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>double</type>
      <name>T</name>
      <anchorfile>classo2scl_1_1rel__boson.html</anchorfile>
      <anchor>aa3b85e7c498e93e1feb728b417f741e1</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>o2scl::rel_fermion</name>
    <filename>classo2scl_1_1rel__fermion.html</filename>
    <base>o2scl::fermion_eval_thermo</base>
    <member kind="function">
      <type></type>
      <name>rel_fermion</name>
      <anchorfile>classo2scl_1_1rel__fermion.html</anchorfile>
      <anchor>acd1f414adf25afebd852c46ecb76b47c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>calc_mu</name>
      <anchorfile>classo2scl_1_1rel__fermion.html</anchorfile>
      <anchor>ac67e3c3cf054b8eb1980e17a7450bb33</anchor>
      <arglist>(fermion &amp;f, double temper)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>calc_density</name>
      <anchorfile>classo2scl_1_1rel__fermion.html</anchorfile>
      <anchor>aedcb1fb028a4eca878e23ee91e8b7f16</anchor>
      <arglist>(fermion &amp;f, double temper)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>pair_mu</name>
      <anchorfile>classo2scl_1_1rel__fermion.html</anchorfile>
      <anchor>a0d6319dd3ba8082ed712fcaa6ef0f12f</anchor>
      <arglist>(fermion &amp;f, double temper)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>pair_density</name>
      <anchorfile>classo2scl_1_1rel__fermion.html</anchorfile>
      <anchor>a4657954489f2866513e9e86c501ab57c</anchor>
      <arglist>(fermion &amp;f, double temper)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual int</type>
      <name>nu_from_n</name>
      <anchorfile>classo2scl_1_1rel__fermion.html</anchorfile>
      <anchor>ad34402527d31ac486c182a3209b9c15f</anchor>
      <arglist>(fermion &amp;f, double temper)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual const char *</type>
      <name>type</name>
      <anchorfile>classo2scl_1_1rel__fermion.html</anchorfile>
      <anchor>a6762673e2092b8709f45fb077aae051a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>err_nonconv</name>
      <anchorfile>classo2scl_1_1rel__fermion.html</anchorfile>
      <anchor>a785a69fa806d9b25a79ed8f2812b5085</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>fermion</type>
      <name>unc</name>
      <anchorfile>classo2scl_1_1rel__fermion.html</anchorfile>
      <anchor>ab5febec8557c1d961d1b3e037a1a3689</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>o2_shared_ptr&lt; inte&lt; funct &gt; &gt;::type</type>
      <name>nit</name>
      <anchorfile>classo2scl_1_1rel__fermion.html</anchorfile>
      <anchor>aef92c35f1a13ff647695895efba693ce</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>o2_shared_ptr&lt; inte&lt; funct &gt; &gt;::type</type>
      <name>dit</name>
      <anchorfile>classo2scl_1_1rel__fermion.html</anchorfile>
      <anchor>ad6adb22693ec28e2f67b6e05f696825e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>o2_shared_ptr&lt; root&lt; funct &gt; &gt;::type</type>
      <name>density_root</name>
      <anchorfile>classo2scl_1_1rel__fermion.html</anchorfile>
      <anchor>a018d39704e044c394cb59f0f61a1aa54</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>min_psi</name>
      <anchorfile>classo2scl_1_1rel__fermion.html</anchorfile>
      <anchor>a4a12fc0a258dd28960da4dbca657c3cb</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>deg_limit</name>
      <anchorfile>classo2scl_1_1rel__fermion.html</anchorfile>
      <anchor>a70ae1aa5b5f3e32bd4b49be4b32ac6f7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>exp_limit</name>
      <anchorfile>classo2scl_1_1rel__fermion.html</anchorfile>
      <anchor>ab20144de680bc2ff2fb09ea4b76e1a37</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>upper_limit_fac</name>
      <anchorfile>classo2scl_1_1rel__fermion.html</anchorfile>
      <anchor>a398fddada96034f44952879f0f824ebd</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>deg_entropy_fac</name>
      <anchorfile>classo2scl_1_1rel__fermion.html</anchorfile>
      <anchor>a23215f2be28272e3c7e9038335f37da1</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" protection="protected">
      <type>double</type>
      <name>density_fun</name>
      <anchorfile>classo2scl_1_1rel__fermion.html</anchorfile>
      <anchor>ae8e4c8087b948d5d587c9090b2a01734</anchor>
      <arglist>(double u)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>double</type>
      <name>energy_fun</name>
      <anchorfile>classo2scl_1_1rel__fermion.html</anchorfile>
      <anchor>a84c440d213c7850190d093a1e4ab00dc</anchor>
      <arglist>(double u)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>double</type>
      <name>entropy_fun</name>
      <anchorfile>classo2scl_1_1rel__fermion.html</anchorfile>
      <anchor>a6834d93824765215504dc7da19fe906f</anchor>
      <arglist>(double u)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>double</type>
      <name>deg_density_fun</name>
      <anchorfile>classo2scl_1_1rel__fermion.html</anchorfile>
      <anchor>aeaab56ddd56a3ebcdf2a430f1d4b1af2</anchor>
      <arglist>(double u)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>double</type>
      <name>deg_energy_fun</name>
      <anchorfile>classo2scl_1_1rel__fermion.html</anchorfile>
      <anchor>a44e1ceebecbba2c717e7c625e2d3e5e6</anchor>
      <arglist>(double u)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>double</type>
      <name>deg_entropy_fun</name>
      <anchorfile>classo2scl_1_1rel__fermion.html</anchorfile>
      <anchor>abd49e59453291d528281102dec2e0c2f</anchor>
      <arglist>(double u)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>double</type>
      <name>solve_fun</name>
      <anchorfile>classo2scl_1_1rel__fermion.html</anchorfile>
      <anchor>ac6f2836572701c0e7a3b3e8f8e233747</anchor>
      <arglist>(double x)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>double</type>
      <name>pair_fun</name>
      <anchorfile>classo2scl_1_1rel__fermion.html</anchorfile>
      <anchor>aa667d278970d33a0a1f5509a7aeff412</anchor>
      <arglist>(double x)</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>double</type>
      <name>T</name>
      <anchorfile>classo2scl_1_1rel__fermion.html</anchorfile>
      <anchor>ae3f97d4eb75b5fa462beb9afd8b73c34</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>fermion *</type>
      <name>fp</name>
      <anchorfile>classo2scl_1_1rel__fermion.html</anchorfile>
      <anchor>af29447021a4ea2df6b1a96f463209ea0</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>o2scl::sn_classical</name>
    <filename>classo2scl_1_1sn__classical.html</filename>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>calc_mu</name>
      <anchorfile>classo2scl_1_1sn__classical.html</anchorfile>
      <anchor>af50477533fa47b579f03192ffb1e1fc7</anchor>
      <arglist>(part_deriv &amp;p, double temper)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>calc_density</name>
      <anchorfile>classo2scl_1_1sn__classical.html</anchorfile>
      <anchor>a96a7036575889ddb22c5113bbbd811f0</anchor>
      <arglist>(part_deriv &amp;p, double temper)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual const char *</type>
      <name>type</name>
      <anchorfile>classo2scl_1_1sn__classical.html</anchorfile>
      <anchor>a568600f6d557158df7fd56d3c41f43ab</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>o2scl::sn_fermion</name>
    <filename>classo2scl_1_1sn__fermion.html</filename>
    <base>o2scl::fermion_deriv_thermo</base>
    <member kind="function">
      <type></type>
      <name>sn_fermion</name>
      <anchorfile>classo2scl_1_1sn__fermion.html</anchorfile>
      <anchor>a383ca0796d7ac4412d967eedeeefd9b5</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>calc_mu</name>
      <anchorfile>classo2scl_1_1sn__fermion.html</anchorfile>
      <anchor>ab037aacb8df48b8b8d44c2d22d019e29</anchor>
      <arglist>(fermion_deriv &amp;f, double temper)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>calc_density</name>
      <anchorfile>classo2scl_1_1sn__fermion.html</anchorfile>
      <anchor>a16df52c687b9bb5d77ad22927a41c774</anchor>
      <arglist>(fermion_deriv &amp;f, double temper)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>pair_mu</name>
      <anchorfile>classo2scl_1_1sn__fermion.html</anchorfile>
      <anchor>ad2f623ae902ea36100270815d151b68b</anchor>
      <arglist>(fermion_deriv &amp;f, double temper)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>pair_density</name>
      <anchorfile>classo2scl_1_1sn__fermion.html</anchorfile>
      <anchor>a025a5cf315e9f670c40e2d9ef2a37ed5</anchor>
      <arglist>(fermion_deriv &amp;f, double temper)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>nu_from_n</name>
      <anchorfile>classo2scl_1_1sn__fermion.html</anchorfile>
      <anchor>a5c5a4771f455bde9dbde6475e88cbb42</anchor>
      <arglist>(fermion_deriv &amp;f, double temper)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_inte</name>
      <anchorfile>classo2scl_1_1sn__fermion.html</anchorfile>
      <anchor>ad592d6949144e9edb9ea0c3bb56a7adf</anchor>
      <arglist>(inte&lt; funct &gt; &amp;unit, inte&lt; funct &gt; &amp;udit)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_density_root</name>
      <anchorfile>classo2scl_1_1sn__fermion.html</anchorfile>
      <anchor>aaedefe2ede51d17748ff0d88b22a60e3</anchor>
      <arglist>(root&lt; funct &gt; &amp;rp)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual const char *</type>
      <name>type</name>
      <anchorfile>classo2scl_1_1sn__fermion.html</anchorfile>
      <anchor>a89a6aed1a1ba826c13bc540f732d8c86</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>deriv_calibrate</name>
      <anchorfile>classo2scl_1_1sn__fermion.html</anchorfile>
      <anchor>a483e89b1dc0b1aa4dba2d8406e391e5f</anchor>
      <arglist>(fermion_deriv &amp;f, int verbose, std::string fname=&quot;&quot;)</arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>exp_limit</name>
      <anchorfile>classo2scl_1_1sn__fermion.html</anchorfile>
      <anchor>aa4daca571c07b71fab6283af2c96fd77</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>deg_limit</name>
      <anchorfile>classo2scl_1_1sn__fermion.html</anchorfile>
      <anchor>acb3316914bb8c2a8e366d3016464c044</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>upper_limit_fac</name>
      <anchorfile>classo2scl_1_1sn__fermion.html</anchorfile>
      <anchor>a1206acdb47c570f3f8de3c68b7a0ff0d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>fermion_deriv</type>
      <name>unc</name>
      <anchorfile>classo2scl_1_1sn__fermion.html</anchorfile>
      <anchor>aeae80a6f3db26add749f9327187297c9</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>inte_qagiu_gsl&lt; funct &gt;</type>
      <name>def_nit</name>
      <anchorfile>classo2scl_1_1sn__fermion.html</anchorfile>
      <anchor>a93610266988cd6440edda1bafcbb60ef</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>inte_qag_gsl&lt; funct &gt;</type>
      <name>def_dit</name>
      <anchorfile>classo2scl_1_1sn__fermion.html</anchorfile>
      <anchor>a8a12d2dee077803f2ae8e7bb5918ddf7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>root_cern&lt; funct &gt;</type>
      <name>def_density_root</name>
      <anchorfile>classo2scl_1_1sn__fermion.html</anchorfile>
      <anchor>a16f471b3585714cc7210291f61100568</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" protection="protected">
      <type>double</type>
      <name>solve_fun</name>
      <anchorfile>classo2scl_1_1sn__fermion.html</anchorfile>
      <anchor>a228012fa0d9edbfa19067ac48aec46a7</anchor>
      <arglist>(double x)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>double</type>
      <name>pair_fun</name>
      <anchorfile>classo2scl_1_1sn__fermion.html</anchorfile>
      <anchor>a3ee448446a000f2ae092c11e36f9c81e</anchor>
      <arglist>(double x)</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>fermion_deriv *</type>
      <name>fp</name>
      <anchorfile>classo2scl_1_1sn__fermion.html</anchorfile>
      <anchor>a14e4778069a30bbc6a57abc0069def12</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>double</type>
      <name>T</name>
      <anchorfile>classo2scl_1_1sn__fermion.html</anchorfile>
      <anchor>ab01b803dd9b3036bac2ac5f43b7a3594</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>inte&lt; funct &gt; *</type>
      <name>nit</name>
      <anchorfile>classo2scl_1_1sn__fermion.html</anchorfile>
      <anchor>ab613cd6bd6040525a2432bd37ad26a19</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>inte&lt; funct &gt; *</type>
      <name>dit</name>
      <anchorfile>classo2scl_1_1sn__fermion.html</anchorfile>
      <anchor>aa584af211e647095ad1d4f9d27c7d8f0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>root&lt; funct &gt; *</type>
      <name>density_root</name>
      <anchorfile>classo2scl_1_1sn__fermion.html</anchorfile>
      <anchor>a5ce70092bbaae893f5c5ce492aa87a9f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>method</name>
      <anchorfile>classo2scl_1_1sn__fermion.html</anchorfile>
      <anchor>a534b1c9533eb5cf59a7e12267f396e83</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const int</type>
      <name>direct</name>
      <anchorfile>classo2scl_1_1sn__fermion.html</anchorfile>
      <anchor>a3be510c41b37d36a122d8fd9b000f73b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const int</type>
      <name>byparts</name>
      <anchorfile>classo2scl_1_1sn__fermion.html</anchorfile>
      <anchor>a911d261d4fd15bd7b628922abb11ffa5</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>o2scl::sn_nr_fermion</name>
    <filename>classo2scl_1_1sn__nr__fermion.html</filename>
    <base>o2scl::fermion_deriv_thermo</base>
    <member kind="function">
      <type></type>
      <name>sn_nr_fermion</name>
      <anchorfile>classo2scl_1_1sn__nr__fermion.html</anchorfile>
      <anchor>a585595c0173549ddae04e364aa7cf0ac</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>calc_mu</name>
      <anchorfile>classo2scl_1_1sn__nr__fermion.html</anchorfile>
      <anchor>ab800aa77b8faf4cbb2ce77d8b4d102ae</anchor>
      <arglist>(fermion_deriv &amp;f, double temper)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>calc_density</name>
      <anchorfile>classo2scl_1_1sn__nr__fermion.html</anchorfile>
      <anchor>a22fc20c7259191fd9e0063d054326696</anchor>
      <arglist>(fermion_deriv &amp;f, double temper)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>pair_mu</name>
      <anchorfile>classo2scl_1_1sn__nr__fermion.html</anchorfile>
      <anchor>a377f77968c6054d776457f8776b1bad7</anchor>
      <arglist>(fermion_deriv &amp;f, double temper)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>pair_density</name>
      <anchorfile>classo2scl_1_1sn__nr__fermion.html</anchorfile>
      <anchor>a5588aec7bdf395a10fe06c77d61d06ae</anchor>
      <arglist>(fermion_deriv &amp;f, double temper)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>nu_from_n</name>
      <anchorfile>classo2scl_1_1sn__nr__fermion.html</anchorfile>
      <anchor>a58e054d8c53e7d9b5d3ec8305e722dbc</anchor>
      <arglist>(fermion_deriv &amp;f, double temper)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_density_root</name>
      <anchorfile>classo2scl_1_1sn__nr__fermion.html</anchorfile>
      <anchor>abbcc53d4058e4e1e0d050c4c8abf6ea9</anchor>
      <arglist>(root&lt; funct &gt; &amp;rp)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual const char *</type>
      <name>type</name>
      <anchorfile>classo2scl_1_1sn__nr__fermion.html</anchorfile>
      <anchor>afccf88f4b261402e49fb0a9a6dbfecdd</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>flimit</name>
      <anchorfile>classo2scl_1_1sn__nr__fermion.html</anchorfile>
      <anchor>ae4a578b707198a5eabc56aa2473a43be</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>fermion_deriv</type>
      <name>unc</name>
      <anchorfile>classo2scl_1_1sn__nr__fermion.html</anchorfile>
      <anchor>a35978fae892938cddda95750fd459cf6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>guess_from_nu</name>
      <anchorfile>classo2scl_1_1sn__nr__fermion.html</anchorfile>
      <anchor>a56bbcf35cd9b967b0b1b31ab75e52616</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>root_cern&lt; funct &gt;</type>
      <name>def_density_root</name>
      <anchorfile>classo2scl_1_1sn__nr__fermion.html</anchorfile>
      <anchor>a90dd19bb79d6220e5122475ec59562dc</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" protection="protected">
      <type>double</type>
      <name>solve_fun</name>
      <anchorfile>classo2scl_1_1sn__nr__fermion.html</anchorfile>
      <anchor>a51b63a09f1e0fbe3006c1073eb836768</anchor>
      <arglist>(double x)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>double</type>
      <name>pair_fun</name>
      <anchorfile>classo2scl_1_1sn__nr__fermion.html</anchorfile>
      <anchor>a367205f59dec74ced506f41fc2ba2d5f</anchor>
      <arglist>(double x)</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>double</type>
      <name>T</name>
      <anchorfile>classo2scl_1_1sn__nr__fermion.html</anchorfile>
      <anchor>a6fe5592c2f8742f8864ca22591b24088</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>fermion_deriv *</type>
      <name>fp</name>
      <anchorfile>classo2scl_1_1sn__nr__fermion.html</anchorfile>
      <anchor>aa77d7cf5c570c35bc8eaa2f36028c7f7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>root&lt; funct &gt; *</type>
      <name>density_root</name>
      <anchorfile>classo2scl_1_1sn__nr__fermion.html</anchorfile>
      <anchor>aa5a97153658c82ce2390cad4387c74f4</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>o2scl::ame_mass</name>
    <filename>classo2scl_1_1ame__mass.html</filename>
    <base>o2scl::nuclear_mass_table</base>
    <class kind="struct">o2scl::ame_mass::entry</class>
    <member kind="function">
      <type></type>
      <name>ame_mass</name>
      <anchorfile>classo2scl_1_1ame__mass.html</anchorfile>
      <anchor>af9a9aa4d9f4bf13e35b8f004710d1905</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual const char *</type>
      <name>type</name>
      <anchorfile>classo2scl_1_1ame__mass.html</anchorfile>
      <anchor>a0b47d9ecc61b7e96aa44484d79ae9655</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>is_included</name>
      <anchorfile>classo2scl_1_1ame__mass.html</anchorfile>
      <anchor>a3bbf7d943cb66051bf91570728ffef35</anchor>
      <arglist>(int Z, int N)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>mass_excess</name>
      <anchorfile>classo2scl_1_1ame__mass.html</anchorfile>
      <anchor>aaee65df7ae93bc2727c3d5e00216d079</anchor>
      <arglist>(int Z, int N)</arglist>
    </member>
    <member kind="function">
      <type>entry</type>
      <name>get_ZN</name>
      <anchorfile>classo2scl_1_1ame__mass.html</anchorfile>
      <anchor>a101dc1e1af2e3ccf95a115ade7342a44</anchor>
      <arglist>(int l_Z, int l_N)</arglist>
    </member>
    <member kind="function">
      <type>entry</type>
      <name>get_ZA</name>
      <anchorfile>classo2scl_1_1ame__mass.html</anchorfile>
      <anchor>a8c6475ed3b8490af8a28b0f0cd1b2f69</anchor>
      <arglist>(int l_Z, int l_A)</arglist>
    </member>
    <member kind="function">
      <type>entry</type>
      <name>get_elA</name>
      <anchorfile>classo2scl_1_1ame__mass.html</anchorfile>
      <anchor>a46ce168d4962e6b25a6c3a7424b3a3f9</anchor>
      <arglist>(std::string l_el, int l_A)</arglist>
    </member>
    <member kind="function">
      <type>entry</type>
      <name>get</name>
      <anchorfile>classo2scl_1_1ame__mass.html</anchorfile>
      <anchor>a205dc1df993dbe2d4692fc3e35798dcb</anchor>
      <arglist>(std::string nucleus)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>is_loaded</name>
      <anchorfile>classo2scl_1_1ame__mass.html</anchorfile>
      <anchor>a72e3d448835c4cab8fe0c449f7a09d79</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>get_nentries</name>
      <anchorfile>classo2scl_1_1ame__mass.html</anchorfile>
      <anchor>a939457d60cce889a3cb61d7873bfa97a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>get_reference</name>
      <anchorfile>classo2scl_1_1ame__mass.html</anchorfile>
      <anchor>ae0359ba204e0a5525adcb3e3fc769e81</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>mass_excess_d</name>
      <anchorfile>classo2scl_1_1nuclear__mass__table.html</anchorfile>
      <anchor>adceeecd759baa5004ba4cbba9c5e23dd</anchor>
      <arglist>(double Z, double N)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual int</type>
      <name>get_nucleus</name>
      <anchorfile>classo2scl_1_1nuclear__mass.html</anchorfile>
      <anchor>ae09d6da7ab9f6254692874996952f5ef</anchor>
      <arglist>(int Z, int N, nucleus &amp;n)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>electron_binding</name>
      <anchorfile>classo2scl_1_1nuclear__mass.html</anchorfile>
      <anchor>a82e8559f6fba9e37c5609071f8a86ee7</anchor>
      <arglist>(double Z)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>binding_energy</name>
      <anchorfile>classo2scl_1_1nuclear__mass.html</anchorfile>
      <anchor>a240471b8e98cbcaa75e6c49036a3658a</anchor>
      <arglist>(int Z, int N)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>binding_energy_d</name>
      <anchorfile>classo2scl_1_1nuclear__mass.html</anchorfile>
      <anchor>ac328add9f2e6978ede52ce3ebfd936b2</anchor>
      <arglist>(double Z, double N)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>total_mass</name>
      <anchorfile>classo2scl_1_1nuclear__mass.html</anchorfile>
      <anchor>a22c3033272eef0347bf084ec20075570</anchor>
      <arglist>(int Z, int N)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>total_mass_d</name>
      <anchorfile>classo2scl_1_1nuclear__mass.html</anchorfile>
      <anchor>a4de19e595362dfe1e3033c732317acdb</anchor>
      <arglist>(double Z, double N)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>atomic_mass</name>
      <anchorfile>classo2scl_1_1nuclear__mass.html</anchorfile>
      <anchor>a4b0104d5111e49f78adc3f3e0bde5a59</anchor>
      <arglist>(int Z, int N)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>atomic_mass_d</name>
      <anchorfile>classo2scl_1_1nuclear__mass.html</anchorfile>
      <anchor>acbe7f41c637118015c22e0e18941d08d</anchor>
      <arglist>(double Z, double N)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>parse_elstring</name>
      <anchorfile>classo2scl_1_1nuclear__mass__info.html</anchorfile>
      <anchor>a11d83323558827c2682fdd0ae28bc1b3</anchor>
      <arglist>(std::string ela, int &amp;Z, int &amp;N, int &amp;A)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>eltoZ</name>
      <anchorfile>classo2scl_1_1nuclear__mass__info.html</anchorfile>
      <anchor>ace3ff294d56921be40769e60a4088923</anchor>
      <arglist>(std::string el)</arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>Ztoel</name>
      <anchorfile>classo2scl_1_1nuclear__mass__info.html</anchorfile>
      <anchor>a9c34e1fbfa0313e205ec603f6f137050</anchor>
      <arglist>(size_t Z)</arglist>
    </member>
    <member kind="function">
      <type>std::string</type>
      <name>tostring</name>
      <anchorfile>classo2scl_1_1nuclear__mass__info.html</anchorfile>
      <anchor>a371a51d2720ee9e8e254640a32109c31</anchor>
      <arglist>(size_t Z, size_t N)</arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const int</type>
      <name>measured</name>
      <anchorfile>classo2scl_1_1ame__mass.html</anchorfile>
      <anchor>a9232f247c824aaf51599d18c524275a6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const int</type>
      <name>estimated</name>
      <anchorfile>classo2scl_1_1ame__mass.html</anchorfile>
      <anchor>a2d46989a32c32427480ada64d2d9aa1d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const int</type>
      <name>not_calculable</name>
      <anchorfile>classo2scl_1_1ame__mass.html</anchorfile>
      <anchor>afb14b029dfdd116a699e1b91804b804a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const int</type>
      <name>intl_computed</name>
      <anchorfile>classo2scl_1_1ame__mass.html</anchorfile>
      <anchor>a4173d1aaa17bd116e75a38bd28861b76</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>int</type>
      <name>n</name>
      <anchorfile>classo2scl_1_1ame__mass.html</anchorfile>
      <anchor>a9595296c1f1392fbb07daadcdf86a487</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>std::string</type>
      <name>reference</name>
      <anchorfile>classo2scl_1_1ame__mass.html</anchorfile>
      <anchor>a40579d336d5fa90908aa23be1949f0b1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>entry *</type>
      <name>mass</name>
      <anchorfile>classo2scl_1_1ame__mass.html</anchorfile>
      <anchor>af3ac86f91ba3325052e6f135f90d2d7d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>int</type>
      <name>last</name>
      <anchorfile>classo2scl_1_1ame__mass.html</anchorfile>
      <anchor>a4a7b5cfab6d434e09551b95c1dd8afef</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>m_neut</name>
      <anchorfile>classo2scl_1_1nuclear__mass.html</anchorfile>
      <anchor>aa48a040dae25fde328ddd462d6631343</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>m_prot</name>
      <anchorfile>classo2scl_1_1nuclear__mass.html</anchorfile>
      <anchor>a805378cdc7ca64b40aba4810c4a83cb8</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>m_elec</name>
      <anchorfile>classo2scl_1_1nuclear__mass.html</anchorfile>
      <anchor>a0fd273b05127be330e2beeaea4f11bb3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>m_amu</name>
      <anchorfile>classo2scl_1_1nuclear__mass.html</anchorfile>
      <anchor>a7cc74aa72c8f20a1835764d89acbb8c3</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef" protection="protected">
      <type>std::map&lt; std::string, int, string_comp &gt;::iterator</type>
      <name>table_it</name>
      <anchorfile>classo2scl_1_1nuclear__mass__info.html</anchorfile>
      <anchor>a86a7779ca994bcd5f8e7b4f4f416c187</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected" static="yes">
      <type>static const int</type>
      <name>nelements</name>
      <anchorfile>classo2scl_1_1nuclear__mass__info.html</anchorfile>
      <anchor>a01589476913ed571a322711ef513becc</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>o2scl::ame_mass::entry</name>
    <filename>structo2scl_1_1ame__mass_1_1entry.html</filename>
    <member kind="variable">
      <type>int</type>
      <name>NMZ</name>
      <anchorfile>structo2scl_1_1ame__mass_1_1entry.html</anchorfile>
      <anchor>a013419f199599603966b197da4a73b23</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>N</name>
      <anchorfile>structo2scl_1_1ame__mass_1_1entry.html</anchorfile>
      <anchor>a7bdeada96739f1f52d847b1a044c5704</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>Z</name>
      <anchorfile>structo2scl_1_1ame__mass_1_1entry.html</anchorfile>
      <anchor>ae74d276fc0d58f939763d530a93634c0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>A</name>
      <anchorfile>structo2scl_1_1ame__mass_1_1entry.html</anchorfile>
      <anchor>afe86f1ceb4d907ac432c03e54a500c2b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>char</type>
      <name>el</name>
      <anchorfile>structo2scl_1_1ame__mass_1_1entry.html</anchorfile>
      <anchor>aca19c27d11d180dda1a73573bef5de19</anchor>
      <arglist>[4]</arglist>
    </member>
    <member kind="variable">
      <type>char</type>
      <name>orig</name>
      <anchorfile>structo2scl_1_1ame__mass_1_1entry.html</anchorfile>
      <anchor>a16550b62f92505d0ea93b78493fb7197</anchor>
      <arglist>[5]</arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>mass</name>
      <anchorfile>structo2scl_1_1ame__mass_1_1entry.html</anchorfile>
      <anchor>a7bbd79ef5b9da67ac3ea9c9abc74729e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>dmass</name>
      <anchorfile>structo2scl_1_1ame__mass_1_1entry.html</anchorfile>
      <anchor>a991357c2f7148f737682db25ee86126d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>mass_acc</name>
      <anchorfile>structo2scl_1_1ame__mass_1_1entry.html</anchorfile>
      <anchor>a1bccd78752a693d0b9b05796632523e1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>be</name>
      <anchorfile>structo2scl_1_1ame__mass_1_1entry.html</anchorfile>
      <anchor>ae46013f5d9331f97426cbf95f3c46966</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>dbe</name>
      <anchorfile>structo2scl_1_1ame__mass_1_1entry.html</anchorfile>
      <anchor>a574282f59149cf7301e9daae1e3f5a38</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>be_acc</name>
      <anchorfile>structo2scl_1_1ame__mass_1_1entry.html</anchorfile>
      <anchor>a7feb1ed9e0f085f9a37f8ce0cea85b07</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>beoa</name>
      <anchorfile>structo2scl_1_1ame__mass_1_1entry.html</anchorfile>
      <anchor>a1bd6751e950bdc7b71dc0196dee8a89e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>dbeoa</name>
      <anchorfile>structo2scl_1_1ame__mass_1_1entry.html</anchorfile>
      <anchor>a3aaf8d186f1b48446c97b8b5b48b42a2</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>beoa_acc</name>
      <anchorfile>structo2scl_1_1ame__mass_1_1entry.html</anchorfile>
      <anchor>a6031700e0def089ed9d17b2947033c8e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>char</type>
      <name>bdmode</name>
      <anchorfile>structo2scl_1_1ame__mass_1_1entry.html</anchorfile>
      <anchor>ae915e3035f048894e4b2245ddea64b71</anchor>
      <arglist>[3]</arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>bde</name>
      <anchorfile>structo2scl_1_1ame__mass_1_1entry.html</anchorfile>
      <anchor>a62f9fdff7a71df325661c40e7931cc17</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>dbde</name>
      <anchorfile>structo2scl_1_1ame__mass_1_1entry.html</anchorfile>
      <anchor>adfc3b08c4f2f191d3ad4553f656eabd5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>bde_acc</name>
      <anchorfile>structo2scl_1_1ame__mass_1_1entry.html</anchorfile>
      <anchor>a1e1d9c9c0b3c24bf45d1a37f6c8a67ec</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>A2</name>
      <anchorfile>structo2scl_1_1ame__mass_1_1entry.html</anchorfile>
      <anchor>a0537f4b06e80b17542ad8b0eb4f3fa61</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>amass</name>
      <anchorfile>structo2scl_1_1ame__mass_1_1entry.html</anchorfile>
      <anchor>a1c4a2cdb471181815cafde743073a8d8</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>damass</name>
      <anchorfile>structo2scl_1_1ame__mass_1_1entry.html</anchorfile>
      <anchor>a857ae1f21bd4b559dce58bfc1c082652</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>amass_acc</name>
      <anchorfile>structo2scl_1_1ame__mass_1_1entry.html</anchorfile>
      <anchor>af107f0a523a5320b189a2616a480ed2a</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>o2scl::ame_mass_exp</name>
    <filename>classo2scl_1_1ame__mass__exp.html</filename>
    <base>o2scl::ame_mass</base>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>is_included</name>
      <anchorfile>classo2scl_1_1ame__mass__exp.html</anchorfile>
      <anchor>a800a770322cb4855811af72f96c0d5f2</anchor>
      <arglist>(int Z, int N)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>o2scl::arb_dist</name>
    <filename>classo2scl_1_1arb__dist.html</filename>
    <base>o2scl::nuclear_dist</base>
    <member kind="function">
      <type></type>
      <name>arb_dist</name>
      <anchorfile>classo2scl_1_1arb__dist.html</anchorfile>
      <anchor>aa0526f8805ab82111306ae35e098cd55</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_dist</name>
      <anchorfile>classo2scl_1_1arb__dist.html</anchorfile>
      <anchor>a0989d91df2c0cbf0785c090349c36383</anchor>
      <arglist>(nuclear_mass &amp;nm, size_t sz, int_vec_t &amp;Zvec, int_vec2_t &amp;Nvec)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add_dist</name>
      <anchorfile>classo2scl_1_1arb__dist.html</anchorfile>
      <anchor>a726823ba92912cf33b2e44301df3c3ce</anchor>
      <arglist>(nuclear_dist &amp;nd)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_dist</name>
      <anchorfile>classo2scl_1_1arb__dist.html</anchorfile>
      <anchor>a878b943f52fe5913e88a414ef7801644</anchor>
      <arglist>(nuclear_mass &amp;nm, std::string expr, int maxA=400, bool include_neutron=false)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual iterator</type>
      <name>begin</name>
      <anchorfile>classo2scl_1_1arb__dist.html</anchorfile>
      <anchor>af1495be9561b4286ad38b8da63b6b8b2</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual iterator</type>
      <name>end</name>
      <anchorfile>classo2scl_1_1arb__dist.html</anchorfile>
      <anchor>a0f5a3794312d5a16fedef5a24b43ebff</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual iterator</type>
      <name>index</name>
      <anchorfile>classo2scl_1_1arb__dist.html</anchorfile>
      <anchor>a941d76be1d9b6e73e0a22b0a88296fbd</anchor>
      <arglist>(size_t i)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual size_t</type>
      <name>size</name>
      <anchorfile>classo2scl_1_1arb__dist.html</anchorfile>
      <anchor>a2c779b210b01f05599e0e9d42f84d4ad</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual nucleus *</type>
      <name>next</name>
      <anchorfile>classo2scl_1_1arb__dist.html</anchorfile>
      <anchor>aeaedb58f84bc338bad9786876493d949</anchor>
      <arglist>(nucleus *np)</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>FunctionParser</type>
      <name>fp</name>
      <anchorfile>classo2scl_1_1arb__dist.html</anchorfile>
      <anchor>a161aecd2f674f1032cb3d6f1178f73c4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>nucleus *</type>
      <name>list</name>
      <anchorfile>classo2scl_1_1arb__dist.html</anchorfile>
      <anchor>ae378f4318c0e9aaa80932bcda8c2dd84</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>size_t</type>
      <name>list_size</name>
      <anchorfile>classo2scl_1_1arb__dist.html</anchorfile>
      <anchor>a554c96847098c838e55168304c93ebdd</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>o2scl::dz_mass_table</name>
    <filename>classo2scl_1_1dz__mass__table.html</filename>
    <base>o2scl::nuclear_mass_table</base>
    <member kind="function">
      <type></type>
      <name>dz_mass_table</name>
      <anchorfile>classo2scl_1_1dz__mass__table.html</anchorfile>
      <anchor>a6844600d6ef5d4fac064799f101265ac</anchor>
      <arglist>(std::string model=&quot;96&quot;, bool external=false)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>is_included</name>
      <anchorfile>classo2scl_1_1dz__mass__table.html</anchorfile>
      <anchor>aa2bf124528faf4b181403572ec6ef276</anchor>
      <arglist>(int Z, int N)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>mass_excess</name>
      <anchorfile>classo2scl_1_1dz__mass__table.html</anchorfile>
      <anchor>af7c159a3ba32f021d2612e4d9f771f25</anchor>
      <arglist>(int Z, int N)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>is_loaded</name>
      <anchorfile>classo2scl_1_1dz__mass__table.html</anchorfile>
      <anchor>a5cc6242a252632bd05dad355d9fb787c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual const char *</type>
      <name>type</name>
      <anchorfile>classo2scl_1_1dz__mass__table.html</anchorfile>
      <anchor>a1832722ffd1ac296abefe061b2d64033</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>get_nentries</name>
      <anchorfile>classo2scl_1_1dz__mass__table.html</anchorfile>
      <anchor>a9ae6e4dc96d0b670394c493ecfbe9a67</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>std::string</type>
      <name>reference</name>
      <anchorfile>classo2scl_1_1dz__mass__table.html</anchorfile>
      <anchor>a777e58b0c9927777b3c55a106589774b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>table</type>
      <name>data</name>
      <anchorfile>classo2scl_1_1dz__mass__table.html</anchorfile>
      <anchor>abb46d47d67f5602cf25efbfd77d825ea</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>int</type>
      <name>last</name>
      <anchorfile>classo2scl_1_1dz__mass__table.html</anchorfile>
      <anchor>a1f3354edea15521a2bc6fd364f1e3086</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>int</type>
      <name>n</name>
      <anchorfile>classo2scl_1_1dz__mass__table.html</anchorfile>
      <anchor>ac36d5ee6df4a9ea7e41c1ca0b4032a80</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>o2scl::dz_mass_fit</name>
    <filename>classo2scl_1_1dz__mass__fit.html</filename>
    <base>o2scl::nuclear_mass_fit</base>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>is_included</name>
      <anchorfile>classo2scl_1_1dz__mass__fit.html</anchorfile>
      <anchor>a7aaa410000a844fb1146a0675a278f7c</anchor>
      <arglist>(int Z, int N)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual const char *</type>
      <name>type</name>
      <anchorfile>classo2scl_1_1dz__mass__fit.html</anchorfile>
      <anchor>a08187d55b3ac4e04e1b6d7c0d27d7ccf</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual int</type>
      <name>fit_fun</name>
      <anchorfile>classo2scl_1_1dz__mass__fit.html</anchorfile>
      <anchor>a60f3b3632f8b296df780b0990f1181eb</anchor>
      <arglist>(size_t nv, const ubvector &amp;x)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual int</type>
      <name>guess_fun</name>
      <anchorfile>classo2scl_1_1dz__mass__fit.html</anchorfile>
      <anchor>a68e70b228b084c58c7cba34bde49a514</anchor>
      <arglist>(size_t nv, ubvector &amp;x)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>binding_energy</name>
      <anchorfile>classo2scl_1_1dz__mass__fit.html</anchorfile>
      <anchor>ae4a5ea410452942b0b7eef5c306dc310</anchor>
      <arglist>(int Z, int N)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>binding_energy_d</name>
      <anchorfile>classo2scl_1_1dz__mass__fit.html</anchorfile>
      <anchor>afd803695ef5278ed8dcbc7145735f3c3</anchor>
      <arglist>(double Z, double N)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>mass_excess</name>
      <anchorfile>classo2scl_1_1dz__mass__fit.html</anchorfile>
      <anchor>af3e3b80f7570565830189fe2371280d0</anchor>
      <arglist>(int Z, int N)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>mass_excess_d</name>
      <anchorfile>classo2scl_1_1dz__mass__fit.html</anchorfile>
      <anchor>aa2e4aa3160e73d5a0ec48ded431f160b</anchor>
      <arglist>(double Z, double N)</arglist>
    </member>
    <member kind="variable">
      <type>ubvector</type>
      <name>b</name>
      <anchorfile>classo2scl_1_1dz__mass__fit.html</anchorfile>
      <anchor>ab307c9bbf4aa06474e288dc139660906</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>size_t</type>
      <name>nfit</name>
      <anchorfile>classo2scl_1_1nuclear__mass__fit.html</anchorfile>
      <anchor>aac1bf2c83bbcd630224f9f8badf22b76</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>tensor3</type>
      <name>onp</name>
      <anchorfile>classo2scl_1_1dz__mass__fit.html</anchorfile>
      <anchor>a67f77e9e6f4bbd711047665f863a9ab3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>ubvector</type>
      <name>y</name>
      <anchorfile>classo2scl_1_1dz__mass__fit.html</anchorfile>
      <anchor>abb40f60a2845d4da86944282f4d2d150</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>ubvector</type>
      <name>pp</name>
      <anchorfile>classo2scl_1_1dz__mass__fit.html</anchorfile>
      <anchor>adcdc81eba99ca6aa4e7a8b06d01faff0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>ubvector</type>
      <name>oei</name>
      <anchorfile>classo2scl_1_1dz__mass__fit.html</anchorfile>
      <anchor>a419ca012d00b38b28b47205091de3096</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>ubvector</type>
      <name>dei</name>
      <anchorfile>classo2scl_1_1dz__mass__fit.html</anchorfile>
      <anchor>afafe08028b21860cc13e057628d4c337</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>ubvector</type>
      <name>qx</name>
      <anchorfile>classo2scl_1_1dz__mass__fit.html</anchorfile>
      <anchor>ae25d33e7f9f117541453eb573d48ea21</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>ubvector</type>
      <name>dx</name>
      <anchorfile>classo2scl_1_1dz__mass__fit.html</anchorfile>
      <anchor>a24813a1298123b65df57d5920f1611bf</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>ubvector</type>
      <name>op</name>
      <anchorfile>classo2scl_1_1dz__mass__fit.html</anchorfile>
      <anchor>a89fcd9846c2119b5f7ac99a1f04f1d90</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>ubvector</type>
      <name>os</name>
      <anchorfile>classo2scl_1_1dz__mass__fit.html</anchorfile>
      <anchor>af0636b9e10e9cfdae327e4b857d38414</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>ubvector</type>
      <name>dyda</name>
      <anchorfile>classo2scl_1_1dz__mass__fit.html</anchorfile>
      <anchor>a94078d4613614c40d58863d0b179f506</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>ubvector_int</type>
      <name>n2</name>
      <anchorfile>classo2scl_1_1dz__mass__fit.html</anchorfile>
      <anchor>afa377f556a07eb51961bd7362b50d9ee</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>ubmatrix_int</type>
      <name>noc</name>
      <anchorfile>classo2scl_1_1dz__mass__fit.html</anchorfile>
      <anchor>a9c593c5ff1ac5cb93f08788f8491db34</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>o2scl::dz_mass_fit_33</name>
    <filename>classo2scl_1_1dz__mass__fit__33.html</filename>
    <base>o2scl::nuclear_mass_fit</base>
    <member kind="function" virtualness="virtual">
      <type>virtual const char *</type>
      <name>type</name>
      <anchorfile>classo2scl_1_1dz__mass__fit__33.html</anchorfile>
      <anchor>a4f6d3c111451d4c5c449fc36d754e62f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual int</type>
      <name>fit_fun</name>
      <anchorfile>classo2scl_1_1dz__mass__fit__33.html</anchorfile>
      <anchor>abcb1d1d428f859898be7bd109ab81a82</anchor>
      <arglist>(size_t nv, const ubvector &amp;x)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual int</type>
      <name>guess_fun</name>
      <anchorfile>classo2scl_1_1dz__mass__fit__33.html</anchorfile>
      <anchor>a3307894869adbd2495f7d1165ada7a17</anchor>
      <arglist>(size_t nv, ubvector &amp;x)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>is_included</name>
      <anchorfile>classo2scl_1_1dz__mass__fit__33.html</anchorfile>
      <anchor>ab2cdafaab9f70adcaf2b6fc3336af47f</anchor>
      <arglist>(int Z, int N)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>binding_energy</name>
      <anchorfile>classo2scl_1_1dz__mass__fit__33.html</anchorfile>
      <anchor>aa11f62b9a4c62ad16789420d4b2f1870</anchor>
      <arglist>(int Z, int N)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>binding_energy_d</name>
      <anchorfile>classo2scl_1_1dz__mass__fit__33.html</anchorfile>
      <anchor>a9410aedeaae6f79f4765d218fbbd8120</anchor>
      <arglist>(double Z, double N)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>mass_excess</name>
      <anchorfile>classo2scl_1_1dz__mass__fit__33.html</anchorfile>
      <anchor>ad460780761e508f838c90b045f3bc884</anchor>
      <arglist>(int Z, int N)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>mass_excess_d</name>
      <anchorfile>classo2scl_1_1dz__mass__fit__33.html</anchorfile>
      <anchor>abc3fa7d8e8a41c11063131bda73e17b0</anchor>
      <arglist>(double Z, double N)</arglist>
    </member>
    <member kind="variable">
      <type>ubvector</type>
      <name>a</name>
      <anchorfile>classo2scl_1_1dz__mass__fit__33.html</anchorfile>
      <anchor>a7654d5219f5680a8070a56349b66f2e9</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>o2scl::frdm_mass</name>
    <filename>classo2scl_1_1frdm__mass.html</filename>
    <base>o2scl::nuclear_mass_fit</base>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>mass_excess_d</name>
      <anchorfile>classo2scl_1_1frdm__mass.html</anchorfile>
      <anchor>aa9512a9346654b2e725b9884b53b571f</anchor>
      <arglist>(double Z, double N)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>mass_excess</name>
      <anchorfile>classo2scl_1_1frdm__mass.html</anchorfile>
      <anchor>ab3cb6d471cf93fa4c52585828c8e3a16</anchor>
      <arglist>(int Z, int N)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual int</type>
      <name>fit_fun</name>
      <anchorfile>classo2scl_1_1frdm__mass.html</anchorfile>
      <anchor>a7f7798337aafac9bf93bb95f0d1e5db9</anchor>
      <arglist>(size_t nv, const ubvector &amp;x)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual int</type>
      <name>guess_fun</name>
      <anchorfile>classo2scl_1_1frdm__mass.html</anchorfile>
      <anchor>a42cfa20b2dfac67a89b488bafa9bad56</anchor>
      <arglist>(size_t nv, ubvector &amp;x)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>drip_binding_energy_d</name>
      <anchorfile>classo2scl_1_1frdm__mass.html</anchorfile>
      <anchor>acc716e4db910d8dc41c0f9611afbd490</anchor>
      <arglist>(double Z, double N, double npout, double nnout, double chi)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>drip_mass_excess_d</name>
      <anchorfile>classo2scl_1_1frdm__mass.html</anchorfile>
      <anchor>ab01fca4c2cf2f9007fd2cc58b21faacd</anchor>
      <arglist>(double Z, double N, double np_out, double nn_out, double chi)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual const char *</type>
      <name>type</name>
      <anchorfile>classo2scl_1_1nuclear__mass__fit.html</anchorfile>
      <anchor>ac2e35c5a908a5acd9d8d71075767e5b1</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>is_included</name>
      <anchorfile>classo2scl_1_1nuclear__mass.html</anchorfile>
      <anchor>a05574bb8df27f7f19de6388c72631d00</anchor>
      <arglist>(int Z, int N)</arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>a1</name>
      <anchorfile>classo2scl_1_1frdm__mass.html</anchorfile>
      <anchor>aa0893508121633d0fb6b8d62c384b5f7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>J</name>
      <anchorfile>classo2scl_1_1frdm__mass.html</anchorfile>
      <anchor>ab91c435470f242f055b19b07c2a57647</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>K</name>
      <anchorfile>classo2scl_1_1frdm__mass.html</anchorfile>
      <anchor>a88de7944ddaef7b2c2c8b4cff2cf7105</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>a2</name>
      <anchorfile>classo2scl_1_1frdm__mass.html</anchorfile>
      <anchor>ac779c6626e7aac8d10771a9af0ea0f0f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>Q</name>
      <anchorfile>classo2scl_1_1frdm__mass.html</anchorfile>
      <anchor>a7a7b59b6dbe3ced048d082ff785f8971</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>a3</name>
      <anchorfile>classo2scl_1_1frdm__mass.html</anchorfile>
      <anchor>a100e803558179f109bef071e8b794bf9</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>ca</name>
      <anchorfile>classo2scl_1_1frdm__mass.html</anchorfile>
      <anchor>af309a5312856f9b6796e9ed9fd229ebc</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>W</name>
      <anchorfile>classo2scl_1_1frdm__mass.html</anchorfile>
      <anchor>a9e9dee7a09712b03799b6d0e254a8fd8</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>ael</name>
      <anchorfile>classo2scl_1_1frdm__mass.html</anchorfile>
      <anchor>af11f1566e05e9b5f7f43f2b8a0784edc</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>rp</name>
      <anchorfile>classo2scl_1_1frdm__mass.html</anchorfile>
      <anchor>ac6c05012ee846e78efd00cf70d48b61c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>r0</name>
      <anchorfile>classo2scl_1_1frdm__mass.html</anchorfile>
      <anchor>a36f603cbf37ce170f18a69cd64dcb6da</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>MH</name>
      <anchorfile>classo2scl_1_1frdm__mass.html</anchorfile>
      <anchor>a7010e27a5e3484fe5ef71b42aee682bd</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>Mn</name>
      <anchorfile>classo2scl_1_1frdm__mass.html</anchorfile>
      <anchor>a928f3a3cf07404f16377030dcd02289f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>e2</name>
      <anchorfile>classo2scl_1_1frdm__mass.html</anchorfile>
      <anchor>a700dfcc5898655d8da2f3d3ffee4d536</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>a</name>
      <anchorfile>classo2scl_1_1frdm__mass.html</anchorfile>
      <anchor>a543a502c1bb18b6d3bd621d6b57076fa</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>aden</name>
      <anchorfile>classo2scl_1_1frdm__mass.html</anchorfile>
      <anchor>a1009d20669a3ed931f82215cb01ecda2</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>rmac</name>
      <anchorfile>classo2scl_1_1frdm__mass.html</anchorfile>
      <anchor>a1e87c1c052ae827216aaec0d98f588ae</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>h</name>
      <anchorfile>classo2scl_1_1frdm__mass.html</anchorfile>
      <anchor>a0ba9455b217661352dd31c79f8086102</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>L</name>
      <anchorfile>classo2scl_1_1frdm__mass.html</anchorfile>
      <anchor>a9d618d7890b9ffab2d0b562121ecdc9b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>C</name>
      <anchorfile>classo2scl_1_1frdm__mass.html</anchorfile>
      <anchor>a09d7d17dbb8a4e6f4c49b98acf0b5eb8</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>gamma</name>
      <anchorfile>classo2scl_1_1frdm__mass.html</anchorfile>
      <anchor>acfb358714b1c5f57a4014a262bb6522e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>amu</name>
      <anchorfile>classo2scl_1_1frdm__mass.html</anchorfile>
      <anchor>a138daed96a5c62fba2602bffe22b1f18</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>nn</name>
      <anchorfile>classo2scl_1_1frdm__mass.html</anchorfile>
      <anchor>adb283fe1bc07e87d898613317b16d279</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>np</name>
      <anchorfile>classo2scl_1_1frdm__mass.html</anchorfile>
      <anchor>a1d186914ff07731cda9db37c2d8e1335</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>Rn</name>
      <anchorfile>classo2scl_1_1frdm__mass.html</anchorfile>
      <anchor>a1400bcbdfe6b442ccebf5947379d1dc2</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>Rp</name>
      <anchorfile>classo2scl_1_1frdm__mass.html</anchorfile>
      <anchor>a441d4e75b3661c86e54fe90d5a261dc4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>double</type>
      <name>kg_to_invfm</name>
      <anchorfile>classo2scl_1_1frdm__mass.html</anchorfile>
      <anchor>a7fc60740b697525157feb3082d918f4f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>double</type>
      <name>Deltap</name>
      <anchorfile>classo2scl_1_1frdm__mass.html</anchorfile>
      <anchor>a30b345536dce8424fbee6317229c15b0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>double</type>
      <name>Deltan</name>
      <anchorfile>classo2scl_1_1frdm__mass.html</anchorfile>
      <anchor>a69f7e53d251da8bf00025d8c5e56d417</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>double</type>
      <name>deltanp</name>
      <anchorfile>classo2scl_1_1frdm__mass.html</anchorfile>
      <anchor>a1c11e36b51bdcaed596f9d03f33c183f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>double</type>
      <name>deltabar</name>
      <anchorfile>classo2scl_1_1frdm__mass.html</anchorfile>
      <anchor>a788aa5c0a3dc70509fd4af66b9fd318f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>double</type>
      <name>epsbar</name>
      <anchorfile>classo2scl_1_1frdm__mass.html</anchorfile>
      <anchor>aa1fc02e8086e5d9079a3b5d6343b6ff2</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>double</type>
      <name>Bs</name>
      <anchorfile>classo2scl_1_1frdm__mass.html</anchorfile>
      <anchor>af428412a1ec4fa936e2a13b179ec9337</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>double</type>
      <name>Bk</name>
      <anchorfile>classo2scl_1_1frdm__mass.html</anchorfile>
      <anchor>abfc7e5bd86e2a1bf574f6390a669efae</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>double</type>
      <name>Br</name>
      <anchorfile>classo2scl_1_1frdm__mass.html</anchorfile>
      <anchor>a047a414035884711e6c1b11af247f1cc</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>double</type>
      <name>Bw</name>
      <anchorfile>classo2scl_1_1frdm__mass.html</anchorfile>
      <anchor>a54a5a6b1c3f6660ea7dc5b48e1e6c3bc</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>double</type>
      <name>Bv</name>
      <anchorfile>classo2scl_1_1frdm__mass.html</anchorfile>
      <anchor>a213e14425c52274be2f86eacc373bcae</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>double</type>
      <name>c1</name>
      <anchorfile>classo2scl_1_1frdm__mass.html</anchorfile>
      <anchor>a2699a5719737e1498850c4ce8bdb3465</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>double</type>
      <name>c2</name>
      <anchorfile>classo2scl_1_1frdm__mass.html</anchorfile>
      <anchor>af4533977083a12d88b2d6fe1e1c18b7b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>double</type>
      <name>c4</name>
      <anchorfile>classo2scl_1_1frdm__mass.html</anchorfile>
      <anchor>a124848612da11a0dea1b32ca9b05bd7d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>double</type>
      <name>c5</name>
      <anchorfile>classo2scl_1_1frdm__mass.html</anchorfile>
      <anchor>a3f1cb25cf19e6a2ba4f0bb9e75329068</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>double</type>
      <name>f0</name>
      <anchorfile>classo2scl_1_1frdm__mass.html</anchorfile>
      <anchor>adc20d76bdca2b8dbed1f80bf34c51eef</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>double</type>
      <name>a0</name>
      <anchorfile>classo2scl_1_1frdm__mass.html</anchorfile>
      <anchor>a69479c17fbe2ab1e46aca0789197631d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>double</type>
      <name>B1</name>
      <anchorfile>classo2scl_1_1frdm__mass.html</anchorfile>
      <anchor>a549e3930db9ef3ebcf9eee694ca4fa53</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>double</type>
      <name>B2</name>
      <anchorfile>classo2scl_1_1frdm__mass.html</anchorfile>
      <anchor>ab48a27484a2d1a7195d7eadc4cfd7e8c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>double</type>
      <name>B3</name>
      <anchorfile>classo2scl_1_1frdm__mass.html</anchorfile>
      <anchor>a4cac487da363e04122df008a4301e9fb</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>double</type>
      <name>B4</name>
      <anchorfile>classo2scl_1_1frdm__mass.html</anchorfile>
      <anchor>a5045c8a2af76fdbcae6d5268b07bc9a6</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>o2scl::mnmsk_mass_entry</name>
    <filename>structo2scl_1_1mnmsk__mass__entry.html</filename>
    <member kind="variable">
      <type>int</type>
      <name>N</name>
      <anchorfile>structo2scl_1_1mnmsk__mass__entry.html</anchorfile>
      <anchor>a035822ef3d2f33f5ddfb602d350c1284</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>Z</name>
      <anchorfile>structo2scl_1_1mnmsk__mass__entry.html</anchorfile>
      <anchor>a41ddee4adccf39b058bbb5f3f05c7094</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>A</name>
      <anchorfile>structo2scl_1_1mnmsk__mass__entry.html</anchorfile>
      <anchor>af69ad6aa131e6e2ee73f763d603a775e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>Emic</name>
      <anchorfile>structo2scl_1_1mnmsk__mass__entry.html</anchorfile>
      <anchor>ae779fa99e64c3efd7fa529b38e042046</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>Mth</name>
      <anchorfile>structo2scl_1_1mnmsk__mass__entry.html</anchorfile>
      <anchor>a1a098b98a60401dd76a3d5f5d740a2be</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>Mexp</name>
      <anchorfile>structo2scl_1_1mnmsk__mass__entry.html</anchorfile>
      <anchor>a78c347c0e234068c0822eeb773b25468</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>sigmaexp</name>
      <anchorfile>structo2scl_1_1mnmsk__mass__entry.html</anchorfile>
      <anchor>a4b30a201f99e55e6971a612ced00e0e3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>EmicFL</name>
      <anchorfile>structo2scl_1_1mnmsk__mass__entry.html</anchorfile>
      <anchor>a2c4608fe0969dbeb3111da5d6e05faf8</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>MthFL</name>
      <anchorfile>structo2scl_1_1mnmsk__mass__entry.html</anchorfile>
      <anchor>abd60c2442f16897dcd5929e2a654374f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>char</type>
      <name>spinp</name>
      <anchorfile>structo2scl_1_1mnmsk__mass__entry.html</anchorfile>
      <anchor>adcb079874162657353876cba09ce6a26</anchor>
      <arglist>[6]</arglist>
    </member>
    <member kind="variable">
      <type>char</type>
      <name>spinn</name>
      <anchorfile>structo2scl_1_1mnmsk__mass__entry.html</anchorfile>
      <anchor>a8fac5989a8619813868919c7488d448c</anchor>
      <arglist>[6]</arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>gapp</name>
      <anchorfile>structo2scl_1_1mnmsk__mass__entry.html</anchorfile>
      <anchor>a333183b3f5576f775f07965b8cb9db6d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>gapn</name>
      <anchorfile>structo2scl_1_1mnmsk__mass__entry.html</anchorfile>
      <anchor>a4a5d48b1977c1a9c35323c4e50c836fe</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>be</name>
      <anchorfile>structo2scl_1_1mnmsk__mass__entry.html</anchorfile>
      <anchor>a5887c8d57f1065de02c670358443d92b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>S1n</name>
      <anchorfile>structo2scl_1_1mnmsk__mass__entry.html</anchorfile>
      <anchor>a9a56933153bc01c5a4021f771d690f7e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>S2n</name>
      <anchorfile>structo2scl_1_1mnmsk__mass__entry.html</anchorfile>
      <anchor>ab7b56671b24b9bc9a9fb7b035b014ef5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>PA</name>
      <anchorfile>structo2scl_1_1mnmsk__mass__entry.html</anchorfile>
      <anchor>a08a5277bd796dea30f935dbaebaf450a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>PAm1</name>
      <anchorfile>structo2scl_1_1mnmsk__mass__entry.html</anchorfile>
      <anchor>ace396fe915c61bdfe6bd284ca9091878</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>PAm2</name>
      <anchorfile>structo2scl_1_1mnmsk__mass__entry.html</anchorfile>
      <anchor>ab0695d07cc417901f1b13936e16bb6a7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>Qbeta</name>
      <anchorfile>structo2scl_1_1mnmsk__mass__entry.html</anchorfile>
      <anchor>a57270cc0398b15ad47c6aabf70e1c04c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>Tbeta</name>
      <anchorfile>structo2scl_1_1mnmsk__mass__entry.html</anchorfile>
      <anchor>a55dcf4661625917bd95952252872f776</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>S1p</name>
      <anchorfile>structo2scl_1_1mnmsk__mass__entry.html</anchorfile>
      <anchor>aced7a8fb519350cd5dd739333c7d3552</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>S2p</name>
      <anchorfile>structo2scl_1_1mnmsk__mass__entry.html</anchorfile>
      <anchor>a583861fa5dfd04b66637a6b0b613daf3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>Qalpha</name>
      <anchorfile>structo2scl_1_1mnmsk__mass__entry.html</anchorfile>
      <anchor>ae7f7984ed0a2c7f21d0d6b3b2277f225</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>Talpha</name>
      <anchorfile>structo2scl_1_1mnmsk__mass__entry.html</anchorfile>
      <anchor>a1f99e4c275cbdb807882f1770bbf5063</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>eps2</name>
      <anchorfile>structo2scl_1_1mnmsk__mass__entry.html</anchorfile>
      <anchor>a603c266b6c378904a35e371b491f8158</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>eps3</name>
      <anchorfile>structo2scl_1_1mnmsk__mass__entry.html</anchorfile>
      <anchor>a3f7ecadc0db7acf6f828a10dfc318131</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>eps4</name>
      <anchorfile>structo2scl_1_1mnmsk__mass__entry.html</anchorfile>
      <anchor>a68469a0d6b7763fb8b0541e5cc53e51d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>eps6</name>
      <anchorfile>structo2scl_1_1mnmsk__mass__entry.html</anchorfile>
      <anchor>ab1877129cd7cfa9ed4f42c73fa5d2350</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>eps6sym</name>
      <anchorfile>structo2scl_1_1mnmsk__mass__entry.html</anchorfile>
      <anchor>a2da872c5530f5d968e5fc2e134997f13</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>beta2</name>
      <anchorfile>structo2scl_1_1mnmsk__mass__entry.html</anchorfile>
      <anchor>a031d96127d71e8a8c4f0829e4f6b3e07</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>beta3</name>
      <anchorfile>structo2scl_1_1mnmsk__mass__entry.html</anchorfile>
      <anchor>afcea246ef117b422ff253f34fcce7b4d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>beta4</name>
      <anchorfile>structo2scl_1_1mnmsk__mass__entry.html</anchorfile>
      <anchor>aedf941d9920f3169f5a1ece5f53a67b7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>beta6</name>
      <anchorfile>structo2scl_1_1mnmsk__mass__entry.html</anchorfile>
      <anchor>aea0dbc43d87d6eb7feeca9a70aea572c</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>o2scl::mnmsk_mass</name>
    <filename>classo2scl_1_1mnmsk__mass.html</filename>
    <base>o2scl::nuclear_mass_table</base>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>is_included</name>
      <anchorfile>classo2scl_1_1mnmsk__mass.html</anchorfile>
      <anchor>ad74d8a7539677f8b6ee0ba26456403d5</anchor>
      <arglist>(int Z, int N)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>mass_excess</name>
      <anchorfile>classo2scl_1_1mnmsk__mass.html</anchorfile>
      <anchor>adf4b26cdb805c530d4f1b2f47122295e</anchor>
      <arglist>(int Z, int N)</arglist>
    </member>
    <member kind="function">
      <type>mnmsk_mass_entry</type>
      <name>get_ZN</name>
      <anchorfile>classo2scl_1_1mnmsk__mass.html</anchorfile>
      <anchor>a06fd40aaa7a530aa5f6d3bc4a005519d</anchor>
      <arglist>(int l_Z, int l_N)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>blank</name>
      <anchorfile>classo2scl_1_1mnmsk__mass.html</anchorfile>
      <anchor>a3df3a8c04dd6fd77df5677b0d74e8444</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>neither</name>
      <anchorfile>classo2scl_1_1mnmsk__mass.html</anchorfile>
      <anchor>a0698778c20c009ca77bd7928a92194cb</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>beta_stable</name>
      <anchorfile>classo2scl_1_1mnmsk__mass.html</anchorfile>
      <anchor>a101d77ebfe3fd1d19fd59c08c0ac85f7</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>beta_plus_and_minus</name>
      <anchorfile>classo2scl_1_1mnmsk__mass.html</anchorfile>
      <anchor>aec0a2b8ab66895e0e6c117c5af54971b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>greater_100</name>
      <anchorfile>classo2scl_1_1mnmsk__mass.html</anchorfile>
      <anchor>a625b22685c0b8ddb1926c503bed072e2</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>very_large</name>
      <anchorfile>classo2scl_1_1mnmsk__mass.html</anchorfile>
      <anchor>ad35d7fb5d500581566cc98068a3aeb7d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual const char *</type>
      <name>type</name>
      <anchorfile>classo2scl_1_1mnmsk__mass.html</anchorfile>
      <anchor>a6b5a88c838b1f993343eb4748b49de21</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>set_data</name>
      <anchorfile>classo2scl_1_1mnmsk__mass.html</anchorfile>
      <anchor>a9b42dee70933618a724343678d8d2a7d</anchor>
      <arglist>(int n_mass, mnmsk_mass_entry *m, std::string ref)</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>int</type>
      <name>n</name>
      <anchorfile>classo2scl_1_1mnmsk__mass.html</anchorfile>
      <anchor>a61c0797c4edfd34e60f1437c5a95c653</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>std::string</type>
      <name>reference</name>
      <anchorfile>classo2scl_1_1mnmsk__mass.html</anchorfile>
      <anchor>a0b9a404884923345909981ae696b10da</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>mnmsk_mass_entry *</type>
      <name>mass</name>
      <anchorfile>classo2scl_1_1mnmsk__mass.html</anchorfile>
      <anchor>a597422165bc34e5d2a82176ed97b7c91</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>int</type>
      <name>last</name>
      <anchorfile>classo2scl_1_1mnmsk__mass.html</anchorfile>
      <anchor>a5382c4d4d4401ff77246d105ed546600</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>o2scl::mnmsk_mass_exp</name>
    <filename>classo2scl_1_1mnmsk__mass__exp.html</filename>
    <base>o2scl::mnmsk_mass</base>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>is_included</name>
      <anchorfile>classo2scl_1_1mnmsk__mass__exp.html</anchorfile>
      <anchor>ae1763431113460d3002152e6186f3f9b</anchor>
      <arglist>(int Z, int N)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>mass_excess</name>
      <anchorfile>classo2scl_1_1mnmsk__mass__exp.html</anchorfile>
      <anchor>ae9a9963814a0b90d9cc810fe3848becc</anchor>
      <arglist>(int Z, int N)</arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>o2scl::hfb_mass_entry</name>
    <filename>structo2scl_1_1hfb__mass__entry.html</filename>
    <member kind="variable">
      <type>int</type>
      <name>N</name>
      <anchorfile>structo2scl_1_1hfb__mass__entry.html</anchorfile>
      <anchor>a4877ce8fdcc537542c8a7390f0cb8695</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>Z</name>
      <anchorfile>structo2scl_1_1hfb__mass__entry.html</anchorfile>
      <anchor>abc7fac23e7c55854efc4e696eb6fbfda</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>A</name>
      <anchorfile>structo2scl_1_1hfb__mass__entry.html</anchorfile>
      <anchor>a3d2d85e59d769c3db7e20491367aa216</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>bet2</name>
      <anchorfile>structo2scl_1_1hfb__mass__entry.html</anchorfile>
      <anchor>a3b81bc6c0f17412664c104c7c29b7793</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>bet4</name>
      <anchorfile>structo2scl_1_1hfb__mass__entry.html</anchorfile>
      <anchor>a946e364b8b93997d068bbd746dad402e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>Rch</name>
      <anchorfile>structo2scl_1_1hfb__mass__entry.html</anchorfile>
      <anchor>ae139cd5d22ac7449e08201529c5cfe58</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>def_wig</name>
      <anchorfile>structo2scl_1_1hfb__mass__entry.html</anchorfile>
      <anchor>a88070e5ac199a21a689c7b7a337e234e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>Sn</name>
      <anchorfile>structo2scl_1_1hfb__mass__entry.html</anchorfile>
      <anchor>a5bf1d43634755cfbf3f64328d05c5aca</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>Sp</name>
      <anchorfile>structo2scl_1_1hfb__mass__entry.html</anchorfile>
      <anchor>a814fb42d626ddf46e57963c146510a7a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>Qbet</name>
      <anchorfile>structo2scl_1_1hfb__mass__entry.html</anchorfile>
      <anchor>a0903206441fd035cbf6a90a1e45acc85</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>Mcal</name>
      <anchorfile>structo2scl_1_1hfb__mass__entry.html</anchorfile>
      <anchor>a50857334e47608cc9ee26eab136da41f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>Err</name>
      <anchorfile>structo2scl_1_1hfb__mass__entry.html</anchorfile>
      <anchor>a6b49eb96c3bbbe5025c7ca6e491872fa</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>o2scl::hfb_sp_mass_entry</name>
    <filename>structo2scl_1_1hfb__sp__mass__entry.html</filename>
    <member kind="variable">
      <type>int</type>
      <name>N</name>
      <anchorfile>structo2scl_1_1hfb__sp__mass__entry.html</anchorfile>
      <anchor>aa19a43aabaff4a9bd08e3ee33f24634d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>Z</name>
      <anchorfile>structo2scl_1_1hfb__sp__mass__entry.html</anchorfile>
      <anchor>ae0d22525a92e3733d58816d3f2991819</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>A</name>
      <anchorfile>structo2scl_1_1hfb__sp__mass__entry.html</anchorfile>
      <anchor>a7e896e970a99ed7f4b2a40c0c4afa369</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>bet2</name>
      <anchorfile>structo2scl_1_1hfb__sp__mass__entry.html</anchorfile>
      <anchor>a981ce10301c890e76deccad47e0df400</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>bet4</name>
      <anchorfile>structo2scl_1_1hfb__sp__mass__entry.html</anchorfile>
      <anchor>a059a950a1192e6f37b6956cd60459056</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>Rch</name>
      <anchorfile>structo2scl_1_1hfb__sp__mass__entry.html</anchorfile>
      <anchor>ad2ed9f09a3c26d9c5f62eec55533df7a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>def_wig</name>
      <anchorfile>structo2scl_1_1hfb__sp__mass__entry.html</anchorfile>
      <anchor>a772fe07210682a74ef1a505e69b0d26d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>Sn</name>
      <anchorfile>structo2scl_1_1hfb__sp__mass__entry.html</anchorfile>
      <anchor>a2f69b9df96a0af6ae78ee9a744b5515c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>Sp</name>
      <anchorfile>structo2scl_1_1hfb__sp__mass__entry.html</anchorfile>
      <anchor>a73b78722811f6cbd0bad18f9e985da77</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>Qbet</name>
      <anchorfile>structo2scl_1_1hfb__sp__mass__entry.html</anchorfile>
      <anchor>a7baaa6ea72f6ba43d0b622705cb11bca</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>Mcal</name>
      <anchorfile>structo2scl_1_1hfb__sp__mass__entry.html</anchorfile>
      <anchor>a9b91d6963628bad1ccd38488569765e4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>Err</name>
      <anchorfile>structo2scl_1_1hfb__sp__mass__entry.html</anchorfile>
      <anchor>adafabc280c0ef6e83a0282f66abc7cbd</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>Jexp</name>
      <anchorfile>structo2scl_1_1hfb__sp__mass__entry.html</anchorfile>
      <anchor>a7a8339fff8bb200ef208a199e0d3346c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>Jth</name>
      <anchorfile>structo2scl_1_1hfb__sp__mass__entry.html</anchorfile>
      <anchor>a0488debd3a4bf3717bed7a9bece477de</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>Pexp</name>
      <anchorfile>structo2scl_1_1hfb__sp__mass__entry.html</anchorfile>
      <anchor>a6c9264629a2de5340d4b34cdc80688c3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>Pth</name>
      <anchorfile>structo2scl_1_1hfb__sp__mass__entry.html</anchorfile>
      <anchor>a2ba21f2ad5f63616947311a63c607ddd</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>o2scl::hfb_mass</name>
    <filename>classo2scl_1_1hfb__mass.html</filename>
    <base>o2scl::nuclear_mass_table</base>
    <member kind="function">
      <type></type>
      <name>hfb_mass</name>
      <anchorfile>classo2scl_1_1hfb__mass.html</anchorfile>
      <anchor>a606fba137434fb2addbff2df1283b50a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>is_included</name>
      <anchorfile>classo2scl_1_1hfb__mass.html</anchorfile>
      <anchor>a2fbc17a3fe4b5a437bb50a6515d15bad</anchor>
      <arglist>(int Z, int N)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>mass_excess</name>
      <anchorfile>classo2scl_1_1hfb__mass.html</anchorfile>
      <anchor>aeb447d21848c60e193d100d7a5fc25a1</anchor>
      <arglist>(int Z, int N)</arglist>
    </member>
    <member kind="function">
      <type>hfb_mass_entry</type>
      <name>get_ZN</name>
      <anchorfile>classo2scl_1_1hfb__mass.html</anchorfile>
      <anchor>acec0c9996200a68a16ab74c1860bc6b2</anchor>
      <arglist>(int l_Z, int l_N)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>is_loaded</name>
      <anchorfile>classo2scl_1_1hfb__mass.html</anchorfile>
      <anchor>a3b786ea4a595c298d126110a18bb5af4</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>blank</name>
      <anchorfile>classo2scl_1_1hfb__mass.html</anchorfile>
      <anchor>a711401b6881cf13c43cc4e4601030914</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual const char *</type>
      <name>type</name>
      <anchorfile>classo2scl_1_1hfb__mass.html</anchorfile>
      <anchor>a01652d8e3e995763a5da2917b950fa09</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>set_data</name>
      <anchorfile>classo2scl_1_1hfb__mass.html</anchorfile>
      <anchor>af11aa12482ca90b9d5e334664d90bf3e</anchor>
      <arglist>(int n_mass, hfb_mass_entry *m, std::string ref)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>get_nentries</name>
      <anchorfile>classo2scl_1_1hfb__mass.html</anchorfile>
      <anchor>a3305bcb029ca73d5c9cc9a07d6165925</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>int</type>
      <name>n</name>
      <anchorfile>classo2scl_1_1hfb__mass.html</anchorfile>
      <anchor>aa352cf0de729718975acc6d722257e04</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>std::string</type>
      <name>reference</name>
      <anchorfile>classo2scl_1_1hfb__mass.html</anchorfile>
      <anchor>a08689767387470ff38707cf29344aad4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>hfb_mass_entry *</type>
      <name>mass</name>
      <anchorfile>classo2scl_1_1hfb__mass.html</anchorfile>
      <anchor>a1ed147d3e5c6f55f3e1bf86caa32296f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>int</type>
      <name>last</name>
      <anchorfile>classo2scl_1_1hfb__mass.html</anchorfile>
      <anchor>aea5a791428654603644dc36d38fdea5a</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>o2scl::hfb_sp_mass</name>
    <filename>classo2scl_1_1hfb__sp__mass.html</filename>
    <base>o2scl::nuclear_mass_table</base>
    <member kind="function">
      <type></type>
      <name>hfb_sp_mass</name>
      <anchorfile>classo2scl_1_1hfb__sp__mass.html</anchorfile>
      <anchor>aadb0c59ebf3b31ae5143365367d96753</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>is_included</name>
      <anchorfile>classo2scl_1_1hfb__sp__mass.html</anchorfile>
      <anchor>a8a2d82de26e2bd33d509fdecfc3dc766</anchor>
      <arglist>(int Z, int N)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>mass_excess</name>
      <anchorfile>classo2scl_1_1hfb__sp__mass.html</anchorfile>
      <anchor>a04693ac1c2c0ec9453f60c89708ac6ca</anchor>
      <arglist>(int Z, int N)</arglist>
    </member>
    <member kind="function">
      <type>hfb_sp_mass_entry</type>
      <name>get_ZN</name>
      <anchorfile>classo2scl_1_1hfb__sp__mass.html</anchorfile>
      <anchor>a7b07914b0846453b3ed00137cc791f4a</anchor>
      <arglist>(int l_Z, int l_N)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual const char *</type>
      <name>type</name>
      <anchorfile>classo2scl_1_1hfb__sp__mass.html</anchorfile>
      <anchor>a67d8ed591a4ecb7c248963f574794131</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>set_data</name>
      <anchorfile>classo2scl_1_1hfb__sp__mass.html</anchorfile>
      <anchor>a63fe5984dddce88c4b3fa920887f1ae6</anchor>
      <arglist>(int n_mass, hfb_sp_mass_entry *m, std::string ref)</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>hfb_sp_mass_entry *</type>
      <name>mass</name>
      <anchorfile>classo2scl_1_1hfb__sp__mass.html</anchorfile>
      <anchor>a6118e7a5bd362ab62bf761ac3a3f30f8</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>int</type>
      <name>n</name>
      <anchorfile>classo2scl_1_1hfb__sp__mass.html</anchorfile>
      <anchor>afadb811d9e88d467f6eac8a775e54015</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>std::string</type>
      <name>reference</name>
      <anchorfile>classo2scl_1_1hfb__sp__mass.html</anchorfile>
      <anchor>a692d636b730c68111677fbbf8d04f41f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>int</type>
      <name>last</name>
      <anchorfile>classo2scl_1_1hfb__sp__mass.html</anchorfile>
      <anchor>a64e51615b10d0a75deab8290935054d9</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="struct">
    <name>o2scl::ktuy_mass_entry</name>
    <filename>structo2scl_1_1ktuy__mass__entry.html</filename>
    <member kind="variable">
      <type>int</type>
      <name>N</name>
      <anchorfile>structo2scl_1_1ktuy__mass__entry.html</anchorfile>
      <anchor>ae95df9854f01e37c07f588f910668ea6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>Z</name>
      <anchorfile>structo2scl_1_1ktuy__mass__entry.html</anchorfile>
      <anchor>a2e4b9bca735a8c847a368e27b34d8c14</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>A</name>
      <anchorfile>structo2scl_1_1ktuy__mass__entry.html</anchorfile>
      <anchor>a11b3ea05c6fcded5569e5ecaf8ec44a7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>Mcal</name>
      <anchorfile>structo2scl_1_1ktuy__mass__entry.html</anchorfile>
      <anchor>a44889c019652cd2f66aab03e6797d845</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>Esh</name>
      <anchorfile>structo2scl_1_1ktuy__mass__entry.html</anchorfile>
      <anchor>ad51a7c48d768e0757993de83a6604979</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>alpha2</name>
      <anchorfile>structo2scl_1_1ktuy__mass__entry.html</anchorfile>
      <anchor>abdd94666d585fe7b09b9d1a929089ca5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>alpha4</name>
      <anchorfile>structo2scl_1_1ktuy__mass__entry.html</anchorfile>
      <anchor>a019f1f255193c0a6c1ad01e4b7569238</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>alpha6</name>
      <anchorfile>structo2scl_1_1ktuy__mass__entry.html</anchorfile>
      <anchor>ad5c71bd4bfc7479233618fbf17c01d4a</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>o2scl::ktuy_mass</name>
    <filename>classo2scl_1_1ktuy__mass.html</filename>
    <base>o2scl::nuclear_mass_table</base>
    <member kind="function">
      <type></type>
      <name>ktuy_mass</name>
      <anchorfile>classo2scl_1_1ktuy__mass.html</anchorfile>
      <anchor>a8dea45dd02e4563495b2bfc91d00ec07</anchor>
      <arglist>(std::string model=&quot;05&quot;, bool external=false)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>is_included</name>
      <anchorfile>classo2scl_1_1ktuy__mass.html</anchorfile>
      <anchor>ac6e0080dc9cc46e78d7416c507123d52</anchor>
      <arglist>(int Z, int N)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>mass_excess</name>
      <anchorfile>classo2scl_1_1ktuy__mass.html</anchorfile>
      <anchor>a23ce78c0c877b9f4dfba67480cae7dd3</anchor>
      <arglist>(int Z, int N)</arglist>
    </member>
    <member kind="function">
      <type>ktuy_mass_entry</type>
      <name>get_ZN</name>
      <anchorfile>classo2scl_1_1ktuy__mass.html</anchorfile>
      <anchor>a06fc69dc7dd346593f1b47fcf53b7c2a</anchor>
      <arglist>(int l_Z, int l_N)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>is_loaded</name>
      <anchorfile>classo2scl_1_1ktuy__mass.html</anchorfile>
      <anchor>a211bb4b6442ca075964b3ba802f5d9e6</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual const char *</type>
      <name>type</name>
      <anchorfile>classo2scl_1_1ktuy__mass.html</anchorfile>
      <anchor>a1ee3cf96672a3c55d384885fe4bbcb7d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>get_nentries</name>
      <anchorfile>classo2scl_1_1ktuy__mass.html</anchorfile>
      <anchor>abbf47340013c4511f9c35f5c34308e03</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>int</type>
      <name>n</name>
      <anchorfile>classo2scl_1_1ktuy__mass.html</anchorfile>
      <anchor>aee6d1df68c2d20a6cfcec623eb636bec</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>std::string</type>
      <name>reference</name>
      <anchorfile>classo2scl_1_1ktuy__mass.html</anchorfile>
      <anchor>a971379479501e8e5f387699f8229cf03</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>ktuy_mass_entry *</type>
      <name>mass</name>
      <anchorfile>classo2scl_1_1ktuy__mass.html</anchorfile>
      <anchor>a5afaf1b2152496938c2513e66d290509</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>int</type>
      <name>last</name>
      <anchorfile>classo2scl_1_1ktuy__mass.html</anchorfile>
      <anchor>ab7f149e31a1252a3ed8c5c565d454b3b</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>o2scl::mass_fit</name>
    <filename>classo2scl_1_1mass__fit.html</filename>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>fit</name>
      <anchorfile>classo2scl_1_1mass__fit.html</anchorfile>
      <anchor>ad372e2c970c71da2b41b914e0d8848d3</anchor>
      <arglist>(nuclear_mass_fit &amp;n, double &amp;res)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>eval</name>
      <anchorfile>classo2scl_1_1mass__fit.html</anchorfile>
      <anchor>ab7bf0b9e743ce1079475ce9fe8291a08</anchor>
      <arglist>(nuclear_mass &amp;n, double &amp;res)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_mmin</name>
      <anchorfile>classo2scl_1_1mass__fit.html</anchorfile>
      <anchor>a48d1c4ce4f20c802b8d198e1feb23536</anchor>
      <arglist>(mmin_base&lt;&gt; &amp;umm)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_dist</name>
      <anchorfile>classo2scl_1_1mass__fit.html</anchorfile>
      <anchor>ac1044247267c78fd2d8a3f1da3ad335f</anchor>
      <arglist>(nuclear_dist &amp;uexp)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_uncerts</name>
      <anchorfile>classo2scl_1_1mass__fit.html</anchorfile>
      <anchor>a934b655f63d5e5a19f0f4419e63a2fe7</anchor>
      <arglist>(size_t nv, vec_t &amp;u)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_exp_mass</name>
      <anchorfile>classo2scl_1_1mass__fit.html</anchorfile>
      <anchor>a0c411574c1bbf9e5757c5d4d696f6b78</anchor>
      <arglist>(nuclear_mass &amp;nm, int maxA=400, bool include_neutron=false)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>eval_isospin_beta</name>
      <anchorfile>classo2scl_1_1mass__fit.html</anchorfile>
      <anchor>a28c26aa43c01bc6b0fa63c05f9c2078a</anchor>
      <arglist>(nuclear_mass &amp;n, ubvector_int &amp;n_qual, ubvector &amp;qual, int max_iso=20)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>eval_isospin</name>
      <anchorfile>classo2scl_1_1mass__fit.html</anchorfile>
      <anchor>a6b0aed8ec0b5f84ba55ac4e914f11092</anchor>
      <arglist>(nuclear_mass &amp;n, ubvector_int &amp;n_qual, ubvector &amp;qual, int min_iso=-8, int max_iso=60)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>min_fun</name>
      <anchorfile>classo2scl_1_1mass__fit.html</anchorfile>
      <anchor>a94895584d93c8ebe5417c9d3cf27cf43</anchor>
      <arglist>(size_t nv, const ubvector &amp;x)</arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>even_even</name>
      <anchorfile>classo2scl_1_1mass__fit.html</anchorfile>
      <anchor>ad9cbbc7d3fb13f83b002891f693617be</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>minZ</name>
      <anchorfile>classo2scl_1_1mass__fit.html</anchorfile>
      <anchor>ad5ad9777816cc261883fdab74304def8</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>minN</name>
      <anchorfile>classo2scl_1_1mass__fit.html</anchorfile>
      <anchor>a23b1a3fd634196b923155319226eb134</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>mmin_simp2</type>
      <name>def_mmin</name>
      <anchorfile>classo2scl_1_1mass__fit.html</anchorfile>
      <anchor>a0df9a0f30f54866aaed85958e9ad60e8</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>full_dist</type>
      <name>def_dist</name>
      <anchorfile>classo2scl_1_1mass__fit.html</anchorfile>
      <anchor>aca9fb46a735b624e891708e0a06a7109</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>ubvector</type>
      <name>uncs</name>
      <anchorfile>classo2scl_1_1mass__fit.html</anchorfile>
      <anchor>af659691b9a4801a202fdb8ffdd384b11</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>mmin_base *</type>
      <name>mm</name>
      <anchorfile>classo2scl_1_1mass__fit.html</anchorfile>
      <anchor>ae1d30d7e0c8ff06f5800e6406b56491c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>nuclear_mass_fit *</type>
      <name>nmf</name>
      <anchorfile>classo2scl_1_1mass__fit.html</anchorfile>
      <anchor>a0f3101fbdefc10e32cf7929b11194448</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>nuclear_dist *</type>
      <name>exp</name>
      <anchorfile>classo2scl_1_1mass__fit.html</anchorfile>
      <anchor>a36f33b9721a36ca15de41926fbb279c6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>fit_method</name>
      <anchorfile>classo2scl_1_1mass__fit.html</anchorfile>
      <anchor>a6260bec9afd453874ddc422d38aa954e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const int</type>
      <name>rms_mass_excess</name>
      <anchorfile>classo2scl_1_1mass__fit.html</anchorfile>
      <anchor>ac394308d1c8c550f0e82b5339e5e4eed</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const int</type>
      <name>rms_binding_energy</name>
      <anchorfile>classo2scl_1_1mass__fit.html</anchorfile>
      <anchor>a4eff6d334d1bac2cc3274591e9c9a111</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const int</type>
      <name>chi_squared_me</name>
      <anchorfile>classo2scl_1_1mass__fit.html</anchorfile>
      <anchor>a7b53e1788037c694a878b076be8ded3a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const int</type>
      <name>chi_squared_be</name>
      <anchorfile>classo2scl_1_1mass__fit.html</anchorfile>
      <anchor>a6bfe4025d4fac1e932ef20b62f72329f</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>o2scl::nuclear_dist</name>
    <filename>classo2scl_1_1nuclear__dist.html</filename>
    <class kind="class">o2scl::nuclear_dist::iterator</class>
    <member kind="function" virtualness="pure">
      <type>virtual iterator</type>
      <name>begin</name>
      <anchorfile>classo2scl_1_1nuclear__dist.html</anchorfile>
      <anchor>a1efe9e868b0e6819e76f3215ca830726</anchor>
      <arglist>()=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual iterator</type>
      <name>end</name>
      <anchorfile>classo2scl_1_1nuclear__dist.html</anchorfile>
      <anchor>a957073e9c71a2a9d4970db1c2e068152</anchor>
      <arglist>()=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual size_t</type>
      <name>size</name>
      <anchorfile>classo2scl_1_1nuclear__dist.html</anchorfile>
      <anchor>a5e598792bdc218b2011ae339ab87c1f7</anchor>
      <arglist>()=0</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="pure">
      <type>virtual nucleus *</type>
      <name>next</name>
      <anchorfile>classo2scl_1_1nuclear__dist.html</anchorfile>
      <anchor>a2b32366d8bb87429a2ce39a76ff1debe</anchor>
      <arglist>(nucleus *np)=0</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>o2scl::nuclear_dist::iterator</name>
    <filename>classo2scl_1_1nuclear__dist_1_1iterator.html</filename>
    <member kind="function">
      <type></type>
      <name>iterator</name>
      <anchorfile>classo2scl_1_1nuclear__dist_1_1iterator.html</anchorfile>
      <anchor>a0a36d6cc60bb1bd8b70497af8109b294</anchor>
      <arglist>(nuclear_dist *ndpp, nucleus *npp)</arglist>
    </member>
    <member kind="function">
      <type>iterator</type>
      <name>operator++</name>
      <anchorfile>classo2scl_1_1nuclear__dist_1_1iterator.html</anchorfile>
      <anchor>ab8542b464b4f6b7cdfeaf7a5d69482a7</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>iterator</type>
      <name>operator++</name>
      <anchorfile>classo2scl_1_1nuclear__dist_1_1iterator.html</anchorfile>
      <anchor>af26205d34d64df924615ee3a1248e0a5</anchor>
      <arglist>(int unused)</arglist>
    </member>
    <member kind="function">
      <type>nucleus *</type>
      <name>operator-&gt;</name>
      <anchorfile>classo2scl_1_1nuclear__dist_1_1iterator.html</anchorfile>
      <anchor>a9578a7024d51c96ed6bf62e8e3adafbd</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>nucleus &amp;</type>
      <name>operator*</name>
      <anchorfile>classo2scl_1_1nuclear__dist_1_1iterator.html</anchorfile>
      <anchor>aa00a70232acefb9cd4f986b5b325ac41</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>nucleus *</type>
      <name>np</name>
      <anchorfile>classo2scl_1_1nuclear__dist_1_1iterator.html</anchorfile>
      <anchor>a86da33f5e67606fc1107bd8494248690</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>nuclear_dist *</type>
      <name>ndp</name>
      <anchorfile>classo2scl_1_1nuclear__dist_1_1iterator.html</anchorfile>
      <anchor>aff608cbb008d615f2c13aa299b468893</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend">
      <type>friend bool</type>
      <name>operator==</name>
      <anchorfile>classo2scl_1_1nuclear__dist_1_1iterator.html</anchorfile>
      <anchor>aeb019e7e65c53440d3e16d7df9ac5ecb</anchor>
      <arglist>(const nuclear_dist::iterator &amp;i1, const nuclear_dist::iterator &amp;i2)</arglist>
    </member>
    <member kind="friend">
      <type>friend bool</type>
      <name>operator!=</name>
      <anchorfile>classo2scl_1_1nuclear__dist_1_1iterator.html</anchorfile>
      <anchor>af56be77fe77aed15f054b4a750c0b5e6</anchor>
      <arglist>(const nuclear_dist::iterator &amp;i1, const nuclear_dist::iterator &amp;i2)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>o2scl::full_dist</name>
    <filename>classo2scl_1_1full__dist.html</filename>
    <base>o2scl::nuclear_dist</base>
    <member kind="function">
      <type></type>
      <name>full_dist</name>
      <anchorfile>classo2scl_1_1full__dist.html</anchorfile>
      <anchor>aea1743d5a97fe29bc5ff815f4a812572</anchor>
      <arglist>(nuclear_mass &amp;nm, int maxA=400, bool include_neutron=false)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>set_dist</name>
      <anchorfile>classo2scl_1_1full__dist.html</anchorfile>
      <anchor>afd9aecf05e9e06c48a972197cdcd4ac2</anchor>
      <arglist>(nuclear_mass &amp;nm, int maxA=400, bool include_neutron=false)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual iterator</type>
      <name>begin</name>
      <anchorfile>classo2scl_1_1full__dist.html</anchorfile>
      <anchor>a15383ff416c70fa46558a2d1038ac9ac</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual iterator</type>
      <name>end</name>
      <anchorfile>classo2scl_1_1full__dist.html</anchorfile>
      <anchor>a775979cce1bf5b863bdc63240f2543c7</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual size_t</type>
      <name>size</name>
      <anchorfile>classo2scl_1_1full__dist.html</anchorfile>
      <anchor>a42108bf5acbd9a8c51f37964ea3df4de</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="protected" virtualness="virtual">
      <type>virtual nucleus *</type>
      <name>next</name>
      <anchorfile>classo2scl_1_1full__dist.html</anchorfile>
      <anchor>a227ab70092c56a0e65ebe3a646859800</anchor>
      <arglist>(nucleus *np)</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>nucleus *</type>
      <name>list</name>
      <anchorfile>classo2scl_1_1full__dist.html</anchorfile>
      <anchor>a1f1857e46443bda9097407afd138be02</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>size_t</type>
      <name>list_size</name>
      <anchorfile>classo2scl_1_1full__dist.html</anchorfile>
      <anchor>a3f2c8708aefbb88c418b656ff04331d6</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>o2scl::nuclear_mass_info</name>
    <filename>classo2scl_1_1nuclear__mass__info.html</filename>
    <member kind="variable" protection="protected">
      <type>std::map&lt; std::string, int, string_comp &gt;</type>
      <name>element_table</name>
      <anchorfile>classo2scl_1_1nuclear__mass__info.html</anchorfile>
      <anchor>a7d47e935a0068028d55130b293f94791</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>std::string</type>
      <name>element_list</name>
      <anchorfile>classo2scl_1_1nuclear__mass__info.html</anchorfile>
      <anchor>a971c3fc469c18ac58e9aeb84be9c59c5</anchor>
      <arglist>[nelements]</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>o2scl::nuclear_mass</name>
    <filename>classo2scl_1_1nuclear__mass.html</filename>
    <base>o2scl::nuclear_mass_info</base>
    <member kind="function" virtualness="virtual">
      <type>virtual const char *</type>
      <name>type</name>
      <anchorfile>classo2scl_1_1nuclear__mass.html</anchorfile>
      <anchor>a724c29a35f55a1fc14e75b96c2947271</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual double</type>
      <name>mass_excess</name>
      <anchorfile>classo2scl_1_1nuclear__mass.html</anchorfile>
      <anchor>aec0f7c0f0a9c796e73fe6b6f2a63f51f</anchor>
      <arglist>(int Z, int N)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual double</type>
      <name>mass_excess_d</name>
      <anchorfile>classo2scl_1_1nuclear__mass.html</anchorfile>
      <anchor>aab3da5836ddb2f78f0ef5173bd836fe1</anchor>
      <arglist>(double Z, double N)=0</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>o2scl::nuclear_mass_table</name>
    <filename>classo2scl_1_1nuclear__mass__table.html</filename>
    <base>o2scl::nuclear_mass</base>
  </compound>
  <compound kind="class">
    <name>o2scl::nuclear_mass_fit</name>
    <filename>classo2scl_1_1nuclear__mass__fit.html</filename>
    <base>o2scl::nuclear_mass</base>
    <member kind="function" virtualness="pure">
      <type>virtual int</type>
      <name>fit_fun</name>
      <anchorfile>classo2scl_1_1nuclear__mass__fit.html</anchorfile>
      <anchor>a6533eeefbb33ca0bc753d1981721f111</anchor>
      <arglist>(size_t nv, const ubvector &amp;x)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual int</type>
      <name>guess_fun</name>
      <anchorfile>classo2scl_1_1nuclear__mass__fit.html</anchorfile>
      <anchor>a8671eaa09c1aa7e40039a01cd11846d9</anchor>
      <arglist>(size_t nv, ubvector &amp;x)=0</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>o2scl::semi_empirical_mass</name>
    <filename>classo2scl_1_1semi__empirical__mass.html</filename>
    <base>o2scl::nuclear_mass_fit</base>
    <member kind="function" virtualness="virtual">
      <type>virtual const char *</type>
      <name>type</name>
      <anchorfile>classo2scl_1_1semi__empirical__mass.html</anchorfile>
      <anchor>a8efc2b563fed5e56373c5562bcf01db1</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>mass_excess_d</name>
      <anchorfile>classo2scl_1_1semi__empirical__mass.html</anchorfile>
      <anchor>abdeacd1af5d8371d55e0ed138e842a55</anchor>
      <arglist>(double Z, double N)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>mass_excess</name>
      <anchorfile>classo2scl_1_1semi__empirical__mass.html</anchorfile>
      <anchor>a7bb255db30c79c10eaaf3c88cb85fb1c</anchor>
      <arglist>(int Z, int N)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual int</type>
      <name>fit_fun</name>
      <anchorfile>classo2scl_1_1semi__empirical__mass.html</anchorfile>
      <anchor>ae986ee45eae21e074e37d6517e12aed9</anchor>
      <arglist>(size_t nv, const ubvector &amp;x)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual int</type>
      <name>guess_fun</name>
      <anchorfile>classo2scl_1_1semi__empirical__mass.html</anchorfile>
      <anchor>a234244a8c37aa280fcad6b5998bf9bc6</anchor>
      <arglist>(size_t nv, ubvector &amp;x)</arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>B</name>
      <anchorfile>classo2scl_1_1semi__empirical__mass.html</anchorfile>
      <anchor>affeaa35834328266ac92646439bbefb9</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>Sv</name>
      <anchorfile>classo2scl_1_1semi__empirical__mass.html</anchorfile>
      <anchor>ac5a57f95ff54bd9a71acbb308b1a4ebd</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>Ss</name>
      <anchorfile>classo2scl_1_1semi__empirical__mass.html</anchorfile>
      <anchor>a55bab813414fd0df4323689aec34b4a5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>Ec</name>
      <anchorfile>classo2scl_1_1semi__empirical__mass.html</anchorfile>
      <anchor>a9a169e4babab0ac07f9627a3defe769a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>Epair</name>
      <anchorfile>classo2scl_1_1semi__empirical__mass.html</anchorfile>
      <anchor>a0445bd0e5d7e3a1988767e7d5e74e431</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>o2scl::ibm_shell_energy</name>
    <filename>classo2scl_1_1ibm__shell__energy.html</filename>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>shell_energy</name>
      <anchorfile>classo2scl_1_1ibm__shell__energy.html</anchorfile>
      <anchor>a7fd7125e9a41c04c35480c1920b1c05d</anchor>
      <arglist>(int Z, int N)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>shell_energy_interp</name>
      <anchorfile>classo2scl_1_1ibm__shell__energy.html</anchorfile>
      <anchor>a72895d3a6785d0227a9966f9574c14c5</anchor>
      <arglist>(double Z, double N)</arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>shells</name>
      <anchorfile>classo2scl_1_1ibm__shell__energy.html</anchorfile>
      <anchor>a69485e8c1aa434ec765fc6875c9fd159</anchor>
      <arglist>[nshells]</arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>shell</name>
      <anchorfile>classo2scl_1_1ibm__shell__energy.html</anchorfile>
      <anchor>af2d276a5753d8afcc465633521922de8</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" static="yes">
      <type>static const size_t</type>
      <name>nshells</name>
      <anchorfile>classo2scl_1_1ibm__shell__energy.html</anchorfile>
      <anchor>a9f66f364e322bce4ca3c0e1a99e2d73e</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>o2scl::dvi_mass</name>
    <filename>classo2scl_1_1dvi__mass.html</filename>
    <base>o2scl::nuclear_mass_fit</base>
    <base>o2scl::ibm_shell_energy</base>
    <member kind="function" virtualness="virtual">
      <type>virtual const char *</type>
      <name>type</name>
      <anchorfile>classo2scl_1_1dvi__mass.html</anchorfile>
      <anchor>ab386db418a829c45f6e94495120387a7</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>mass_excess_d</name>
      <anchorfile>classo2scl_1_1dvi__mass.html</anchorfile>
      <anchor>a7365e8116cd1b47b02e159ff3795dafb</anchor>
      <arglist>(double Z, double N)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual double</type>
      <name>mass_excess</name>
      <anchorfile>classo2scl_1_1dvi__mass.html</anchorfile>
      <anchor>a78a56d4342abf549d5645b8d39f997a3</anchor>
      <arglist>(int Z, int N)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual int</type>
      <name>fit_fun</name>
      <anchorfile>classo2scl_1_1dvi__mass.html</anchorfile>
      <anchor>a5a60c190e690154aeea081669c516353</anchor>
      <arglist>(size_t nv, const ubvector &amp;x)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual int</type>
      <name>guess_fun</name>
      <anchorfile>classo2scl_1_1dvi__mass.html</anchorfile>
      <anchor>ad23ea54ce4caf07162ab028971ec0014</anchor>
      <arglist>(size_t nv, ubvector &amp;x)</arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>av</name>
      <anchorfile>classo2scl_1_1dvi__mass.html</anchorfile>
      <anchor>a52020fc35fd11926c0f796667665af11</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>as</name>
      <anchorfile>classo2scl_1_1dvi__mass.html</anchorfile>
      <anchor>a682f0d225403e3d884b3f383f467a800</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>sv</name>
      <anchorfile>classo2scl_1_1dvi__mass.html</anchorfile>
      <anchor>acd27896335dea0b92cf2d0e0b1575844</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>ac</name>
      <anchorfile>classo2scl_1_1dvi__mass.html</anchorfile>
      <anchor>a3bd076bc5c8fa55eaace87d4ceab4e80</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>ap</name>
      <anchorfile>classo2scl_1_1dvi__mass.html</anchorfile>
      <anchor>aeaa583178464b2dd59a2462a13472b35</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>y</name>
      <anchorfile>classo2scl_1_1dvi__mass.html</anchorfile>
      <anchor>a660afce4675d580531a4c41858af1de7</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>o2scl::rms_radius</name>
    <filename>classo2scl_1_1rms__radius.html</filename>
    <member kind="function">
      <type>int</type>
      <name>eval_rms_rho</name>
      <anchorfile>classo2scl_1_1rms__radius.html</anchorfile>
      <anchor>a62c086b7fd4551c154d3835847561963</anchor>
      <arglist>(double rho0, double N, double d, double &amp;Rcd, double &amp;Rfermi, double &amp;Rrms)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>eval_rms_rsq</name>
      <anchorfile>classo2scl_1_1rms__radius.html</anchorfile>
      <anchor>a2cc15f31d89e3efb18ae2d034d010d6c</anchor>
      <arglist>(double Rfermi, double N, double d, double &amp;rho0, double &amp;Rcd, double &amp;Rrms)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>double</type>
      <name>iand</name>
      <anchorfile>classo2scl_1_1rms__radius.html</anchorfile>
      <anchor>a59d7745bb57c4f8cb4ac7d74f79d6098</anchor>
      <arglist>(double r)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>double</type>
      <name>iand2</name>
      <anchorfile>classo2scl_1_1rms__radius.html</anchorfile>
      <anchor>a77b2902b92cc7be03da29c67b118483b</anchor>
      <arglist>(double r)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>double</type>
      <name>solve</name>
      <anchorfile>classo2scl_1_1rms__radius.html</anchorfile>
      <anchor>a2358f011a6c02086daaa632096de97c8</anchor>
      <arglist>(double x)</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>double</type>
      <name>urho0</name>
      <anchorfile>classo2scl_1_1rms__radius.html</anchorfile>
      <anchor>aad7f96e4f2953a1e8a0908e17511cfdd</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>double</type>
      <name>ud</name>
      <anchorfile>classo2scl_1_1rms__radius.html</anchorfile>
      <anchor>aac5c1ed6bcac9b9626d91d1be404a0e4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>double</type>
      <name>uRfermi</name>
      <anchorfile>classo2scl_1_1rms__radius.html</anchorfile>
      <anchor>a5ef03e737014de72721a2024081be262</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>double</type>
      <name>uN</name>
      <anchorfile>classo2scl_1_1rms__radius.html</anchorfile>
      <anchor>accaa81e0d2010b19807eff76b2c00e25</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>inte_qagiu_gsl&lt; funct &gt;</type>
      <name>it</name>
      <anchorfile>classo2scl_1_1rms__radius.html</anchorfile>
      <anchor>a4fe7e101f36d984afeaaf9c6607ca2a3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>root_cern&lt; funct &gt;</type>
      <name>cr</name>
      <anchorfile>classo2scl_1_1rms__radius.html</anchorfile>
      <anchor>ac22fed64b21c9819035087729d209da3</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>o2scl::nucleus</name>
    <filename>classo2scl_1_1nucleus.html</filename>
    <base>o2scl::part</base>
    <member kind="variable">
      <type>int</type>
      <name>Z</name>
      <anchorfile>classo2scl_1_1nucleus.html</anchorfile>
      <anchor>a2f9bbfb9b0a4ae4036ed4218d31db1d7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>N</name>
      <anchorfile>classo2scl_1_1nucleus.html</anchorfile>
      <anchor>ac84d317bc503ea685d93ee1f7b74f6ec</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>A</name>
      <anchorfile>classo2scl_1_1nucleus.html</anchorfile>
      <anchor>a1b5323dcdbcbccf8a3cc2431b4caeeb9</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>mex</name>
      <anchorfile>classo2scl_1_1nucleus.html</anchorfile>
      <anchor>a5424d65276cf2a60114d0fa60f3998b2</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>be</name>
      <anchorfile>classo2scl_1_1nucleus.html</anchorfile>
      <anchor>a51a8cc89cb6e43161be1fce2ff4d7f72</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>o2scl::nuclear_reaction</name>
    <filename>classo2scl_1_1nuclear__reaction.html</filename>
    <member kind="function">
      <type>std::string</type>
      <name>to_string</name>
      <anchorfile>classo2scl_1_1nuclear__reaction.html</anchorfile>
      <anchor>ac0c34e123f0a13141d0ad9e90f3716de</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>clear</name>
      <anchorfile>classo2scl_1_1nuclear__reaction.html</anchorfile>
      <anchor>a752da1d766232e13bdb8e710e99990fd</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>nuclear_reaction</name>
      <anchorfile>classo2scl_1_1nuclear__reaction.html</anchorfile>
      <anchor>a90b7772dfa8f83246951b3f3619055af</anchor>
      <arglist>(const nuclear_reaction &amp;nr)</arglist>
    </member>
    <member kind="function">
      <type>nuclear_reaction &amp;</type>
      <name>operator=</name>
      <anchorfile>classo2scl_1_1nuclear__reaction.html</anchorfile>
      <anchor>a11fd287089f57fd8c574c3c259b3dcba</anchor>
      <arglist>(const nuclear_reaction &amp;nr)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>rate</name>
      <anchorfile>classo2scl_1_1nuclear__reaction.html</anchorfile>
      <anchor>a6b9d68374a765a0d9d254fed7e3f7bf0</anchor>
      <arglist>(double T9)</arglist>
    </member>
    <member kind="variable">
      <type>size_t</type>
      <name>chap</name>
      <anchorfile>classo2scl_1_1nuclear__reaction.html</anchorfile>
      <anchor>a009ce31f3450950a97c7519888301a86</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>std::string</type>
      <name>name</name>
      <anchorfile>classo2scl_1_1nuclear__reaction.html</anchorfile>
      <anchor>a1541c27dab57e7837f10642339787bef</anchor>
      <arglist>[6]</arglist>
    </member>
    <member kind="variable">
      <type>std::string</type>
      <name>ref</name>
      <anchorfile>classo2scl_1_1nuclear__reaction.html</anchorfile>
      <anchor>a4778c5e3ff641c5e0108bc87ff9b73e3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>char</type>
      <name>type</name>
      <anchorfile>classo2scl_1_1nuclear__reaction.html</anchorfile>
      <anchor>a62ddb6c5beef9231ea7e2b6586135e38</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>char</type>
      <name>rev</name>
      <anchorfile>classo2scl_1_1nuclear__reaction.html</anchorfile>
      <anchor>a12197c10f36c1635bbe5b45f01557a9d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>Q</name>
      <anchorfile>classo2scl_1_1nuclear__reaction.html</anchorfile>
      <anchor>a43990fa2b900f1ce3b5b3122948562a7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>a</name>
      <anchorfile>classo2scl_1_1nuclear__reaction.html</anchorfile>
      <anchor>a732ca65513fdcd50ff1b2f6e93231acd</anchor>
      <arglist>[7]</arglist>
    </member>
    <member kind="variable">
      <type>size_t</type>
      <name>Z</name>
      <anchorfile>classo2scl_1_1nuclear__reaction.html</anchorfile>
      <anchor>a4cb6f97dbf4235a6cd89e1590b58bded</anchor>
      <arglist>[6]</arglist>
    </member>
    <member kind="variable">
      <type>size_t</type>
      <name>A</name>
      <anchorfile>classo2scl_1_1nuclear__reaction.html</anchorfile>
      <anchor>a0911049115186ff7d2358887cf881623</anchor>
      <arglist>[6]</arglist>
    </member>
    <member kind="variable">
      <type>size_t</type>
      <name>isomer</name>
      <anchorfile>classo2scl_1_1nuclear__reaction.html</anchorfile>
      <anchor>a55472dfa8b2873fcdc1c40953c26b2ae</anchor>
      <arglist>[6]</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>o2scl::reaction_lib</name>
    <filename>classo2scl_1_1reaction__lib.html</filename>
    <member kind="function">
      <type>int</type>
      <name>read_file_reaclib2</name>
      <anchorfile>classo2scl_1_1reaction__lib.html</anchorfile>
      <anchor>afdeb53c4f719cb25ede62de2c0f5a551</anchor>
      <arglist>(std::string fname)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>find_in_chap</name>
      <anchorfile>classo2scl_1_1reaction__lib.html</anchorfile>
      <anchor>abfa07d12fca90dce530a9381ab4fef25</anchor>
      <arglist>(std::vector&lt; nuclear_reaction &gt; &amp;nrl, size_t chap, std::string nuc1, std::string nuc2=&quot;&quot;, std::string nuc3=&quot;&quot;, std::string nuc4=&quot;&quot;, std::string nuc5=&quot;&quot;, std::string nuc6=&quot;&quot;)</arglist>
    </member>
    <member kind="variable">
      <type>std::vector&lt; nuclear_reaction &gt;</type>
      <name>lib</name>
      <anchorfile>classo2scl_1_1reaction__lib.html</anchorfile>
      <anchor>a5dc079a78413703ba64ec7c544680fe3</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" protection="protected">
      <type>bool</type>
      <name>matches</name>
      <anchorfile>classo2scl_1_1reaction__lib.html</anchorfile>
      <anchor>a3e3cf22d313d8088713f553b8cf25618</anchor>
      <arglist>(size_t ul, size_t ri)</arglist>
    </member>
  </compound>
  <compound kind="dir">
    <name>nuclei</name>
    <path>/Users/awsteiner/svn/osf/trunk/src/nuclei/</path>
    <filename>dir_0ea211a22ec1925bcd4ff15e6ee002ab.html</filename>
    <file>ame_mass.h</file>
    <file>arb_dist.h</file>
    <file>dz_mass.h</file>
    <file>frdm_mass.h</file>
    <file>hdf_nucmass_io.h</file>
    <file>hfb_mass.h</file>
    <file>ktuy_mass.h</file>
    <file>mass_fit.h</file>
    <file>nuclear_dist.h</file>
    <file>nuclear_mass.h</file>
    <file>nucleus.h</file>
    <file>reaction_lib.h</file>
  </compound>
  <compound kind="dir">
    <name>part</name>
    <path>/Users/awsteiner/svn/osf/trunk/src/part/</path>
    <filename>dir_6cf79716bc08cf8402a4485dbab34a0b.html</filename>
    <file>boson.h</file>
    <file>classical.h</file>
    <file>eff_boson.h</file>
    <file>eff_fermion.h</file>
    <file>fermion.h</file>
    <file>mag_fermion_zerot.h</file>
    <file>nonrel_fermion.h</file>
    <file>part.h</file>
    <file>part_deriv.h</file>
    <file>quark.h</file>
    <file>rel_boson.h</file>
    <file>rel_fermion.h</file>
    <file>sn_classical.h</file>
    <file>sn_fermion.h</file>
    <file>sn_nr_fermion.h</file>
  </compound>
  <compound kind="dir">
    <name>src</name>
    <path>/Users/awsteiner/svn/osf/trunk/src/</path>
    <filename>dir_68267d1309a1af8e8297ef4c3efbcdba.html</filename>
    <dir>nuclei</dir>
    <dir>part</dir>
  </compound>
  <compound kind="page">
    <name>index</name>
    <title>O2scl Particles and Nuclei Sub-Library User&apos;s Guide</title>
    <filename>index</filename>
    <docanchor file="index" title="Feature Overview">partfeat_section</docanchor>
    <docanchor file="index" title="Quick Reference to User&apos;s Guide">prt_ug_section</docanchor>
  </compound>
</tagfile>
