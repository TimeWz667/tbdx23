const { description } = require('../../package')

module.exports = {
  title: 'Cascade tester',

  description: description,
  publicPath: '/tbdx23/',
  head: [
    ['meta', { name: 'theme-color', content: '#3eaf7c' }],
    ['meta', { name: 'apple-mobile-web-app-capable', content: 'yes' }],
    ['meta', { name: 'apple-mobile-web-app-status-bar-style', content: 'black' }]
  ],

  /**
   * Theme configuration, here is the default theme configuration for VuePress.
   *
   * ref：https://v1.vuepress.vuejs.org/theme/default-theme-config.html
   */
  themeConfig: {
    repo: '',
    editLinks: false,
    docsDir: '',
    editLinkText: '',
    lastUpdated: false,
    nav: [
      {
        text: 'Home',
        link: '/',
      },
      {
        text: 'India',
        link: '/ind/',
      },
      {
        text: 'South Africa',
        link: '/rsa/'
      }
    ],
    sidebar: {
      '/ind/': [
        '',
        'inputs',
      ],
      '/rsa/': [
        '',
        'inputs'
      ],
    },
    sidebarDepth: 3
  },

  /**
   * Apply plugins，ref：https://v1.vuepress.vuejs.org/zh/plugin/
   */
  plugins: [
    '@vuepress/plugin-back-to-top',
    '@vuepress/plugin-medium-zoom',
  ]
}
