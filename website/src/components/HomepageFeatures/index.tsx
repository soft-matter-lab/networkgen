import type {ReactNode} from 'react';
import useBaseUrl from '@docusaurus/useBaseUrl';
import Heading from '@theme/Heading';
import styles from './styles.module.css';

type FeatureItem = {
  title: string;
  imgSrc?: string;
  imgBg?: string;
  imgNode?: ReactNode;
  description: ReactNode;
};

function LammpsImageCard() {
  return (
    <div className={styles.lammpsImg} style={{margin: '-10px -5px 0px -5px'}}>
      <div className={styles.lammpsWordmark}>LAMMPS</div>
      <div className={styles.lammpsSubtext}>Large-scale Atomic/Molecular<br/>Massively Parallel Simulator</div>
      <div className={styles.lammpsChips}>
        <code className={styles.chip}>PolyNetwork.data</code>
        <code className={styles.chip}>bond_table.txt</code>
        <code className={styles.chip}>PolyVisual.lammpstrj</code>
      </div>
    </div>
  );
}

const FeatureList: FeatureItem[] = [
  {
    title: 'Control geometry & architecture',
    description: (
      <>
        Generate random or hexagonal lattice networks with tunable,
        node density, and per-atom bond limits. Full control over domain size
        and boundary conditions.
      </>
    ),
  },
  {
    title: 'Precise network statistics',
    description: (
      <>
        Precisely control strand length distributions — monodisperse,
        polydisperse, or bimodal — with full parameterisation of each mode
        including prestretch and Kuhn length assignment.
      </>
    ),
  },
  {
    title: 'Defect and disorder generation',
    description: (
      <>
        Introduce lattice disorder, create defects, or network heterogeneities. 
        Use Gaussian, fixed, or exponential size distributions.
        Control void placement, shape roughness, clustering, and more.
      </>
    ),
  },
  {
    title: 'Deploy to LAMMPS',
    imgNode: <LammpsImageCard />,
    description: (
      <>
        NetworkGen data and table files are designed to work directly with LAMMPS. 
        Generate your network and get started with your simulations!
      </>
    ),
  },
];

function Feature({title, imgSrc, imgBg, imgNode, description}: FeatureItem) {
  return (
    <div className={styles.featureCard}>
      <div className={styles.featureImgWrap} style={imgBg ? {background: imgBg} : undefined}>
        {imgNode
          ? imgNode
          : imgSrc
            ? <img src={imgSrc} alt={title} className={styles.featureImg} />
            : null
        }
      </div>
      <div className={styles.featureContent}>
        <Heading as="h3" className={styles.featureTitle}>{title}</Heading>
        <p className={styles.featureDesc}>{description}</p>
      </div>
    </div>
  );
}

export default function HomepageFeatures(): ReactNode {
  const archImg = useBaseUrl('/img/network_arch_main.svg');
  const gaussImg = useBaseUrl('/img/network_gauss_main2.png');
  const defectImg = useBaseUrl('/img/defect_network.svg');

  const withImages: FeatureItem[] = FeatureList.map((f, i) => {
    if (i === 0) return {...f, imgSrc: archImg, imgBg: '#1a1628'};
    if (i === 1) return {...f, imgSrc: gaussImg, imgBg: '#1a1628'};
    if (i === 2) return {...f, imgSrc: defectImg, imgBg: '#1a1628'};
    return f;
  });

  return (
    <section className={styles.features}>
      <div className="container">
        <div className={styles.featureGrid}>
          {withImages.map((props, idx) => (
            <Feature key={idx} {...props} />
          ))}
        </div>
      </div>
    </section>
  );
}
