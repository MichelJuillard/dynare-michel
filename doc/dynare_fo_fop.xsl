<?xml version='1.0'?>
<xsl:stylesheet xmlns:xsl="http://www.w3.org/1999/XSL/Transform"
                xmlns:fo="http://www.w3.org/1999/XSL/Format"
                xmlns:doc="http://nwalsh.com/xsl/documentation/1.0"
                exclude-result-prefixes="doc"
                version='1.0'>

<!-- It is important to use indent="no" here, otherwise verbatim -->
<!-- environments get broken by indented tags...at least when the -->
<!-- callout extension is used...at least with some processors -->
<!-- <xsl:import href="c:/cygwin/usr/share/docbook-xsl/fo/docbook.xsl"/> -->
<xsl:import href="e:/docbook-xsl-1.65.1/fo/docbook.xsl"/>

<xsl:param name="use.extensions">1</xsl:param>
<xsl:param name="fop.extensions">1</xsl:param>
<xsl:param name="refentry.generate.name">0</xsl:param>
<xsl:param name="refentry.generate.title">1</xsl:param>
<xsl:template match="refsynopsisdiv">
  <fo:block font-weight="bold" font-size="16pt" font-family="sans-serif">
    <xsl:call-template name="gentext">
                  <xsl:with-param name="key" select="'RefSynopsisDiv'"/>
    </xsl:call-template>
  </fo:block>
    <xsl:apply-templates/>
</xsl:template>

<xsl:template match="refnamediv">
  <xsl:variable name="section.level">
    <xsl:call-template name="refentry.level">
      <xsl:with-param name="node" select="ancestor::refentry"/>
    </xsl:call-template>
  </xsl:variable>

  <xsl:variable name="reftitle">
    <xsl:choose>
      <xsl:when test="$refentry.generate.name != 0">
        <xsl:call-template name="gentext">
          <xsl:with-param name="key" select="'RefName'"/>
        </xsl:call-template>
      </xsl:when>
      <xsl:when test="$refentry.generate.title != 0">
        <xsl:choose>
          <xsl:when test="../refmeta/refentrytitle">
            <xsl:apply-templates select="../refmeta/refentrytitle"/>
          </xsl:when>
          <xsl:otherwise>
            <xsl:apply-templates select="refname[1]"/>
          </xsl:otherwise>
        </xsl:choose>
      </xsl:when>
    </xsl:choose>
  </xsl:variable>

  <!-- xsl:use-attribute-sets takes only a Qname, not a variable -->
  <fo:block>
    <xsl:choose>
      <xsl:when test="$section.level = 1">
        <fo:block xsl:use-attribute-sets="refentry.title.properties">
          <fo:block xsl:use-attribute-sets="section.title.level1.properties">
            <xsl:value-of select="$reftitle"/>
          </fo:block>
        </fo:block>
      </xsl:when>
      <xsl:when test="$section.level = 2">
        <fo:block xsl:use-attribute-sets="refentry.title.properties">
          <fo:block xsl:use-attribute-sets="section.title.level2.properties">
            <xsl:value-of select="$reftitle"/>
          </fo:block>
        </fo:block>
      </xsl:when>
      <xsl:when test="$section.level = 3">
        <fo:block xsl:use-attribute-sets="refentry.title.properties">
          <fo:block xsl:use-attribute-sets="section.title.level3.properties">
            <xsl:value-of select="$reftitle"/>
          </fo:block>
        </fo:block>
      </xsl:when>
      <xsl:when test="$section.level = 4">
        <fo:block xsl:use-attribute-sets="refentry.title.properties">
          <fo:block xsl:use-attribute-sets="section.title.level4.properties">
            <xsl:value-of select="$reftitle"/>
          </fo:block>
        </fo:block>
      </xsl:when>
      <xsl:when test="$section.level = 5">
        <fo:block xsl:use-attribute-sets="refentry.title.properties">
          <fo:block xsl:use-attribute-sets="section.title.level5.properties">
            <xsl:value-of select="$reftitle"/>
          </fo:block>
        </fo:block>
      </xsl:when>
      <xsl:otherwise>
        <fo:block xsl:use-attribute-sets="refentry.title.properties">
          <fo:block xsl:use-attribute-sets="section.title.level6.properties">
            <xsl:value-of select="$reftitle"/>
          </fo:block>
        </fo:block>
      </xsl:otherwise>
    </xsl:choose>

    <fo:block space-after="1em">
      <xsl:choose>
        <xsl:when test="../refmeta/refentrytitle">
          <xsl:apply-templates select="../refmeta/refentrytitle"/>
        </xsl:when>
        <xsl:otherwise>
          <xsl:apply-templates select="refname[1]"/>
        </xsl:otherwise>
      </xsl:choose>
      <xsl:apply-templates select="refpurpose"/>
    </fo:block>
  </fo:block>
  
</xsl:template>
</xsl:stylesheet>
